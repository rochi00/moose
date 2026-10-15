//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVPropertyInterface.h"

#include "KokkosConstantFunction.h"
#include "KokkosLinearFVBoundaryCondition.h"
#include "KokkosLinearFVVariableInterface.h"
#include "KokkosParsedFunction.h"
#include "KokkosParsedFunctorMaterial.h"
#include "KokkosReferenceWrapper.h"

#include "Attributes.h"
#include "FEProblemBase.h"
#include "FVGeometricAverage.h"
#include "FVHarmonicAverage.h"
#include "LinearFVBoundaryCondition.h"
#include "MooseLinearVariableFV.h"
#include "MooseObject.h"
#include "MooseUtils.h"

#include <algorithm>
#include <limits>

namespace
{
using namespace Moose::Kokkos;

/// Interpolations of the cell values of a property to the faces with the variable on both sides
enum class Interpolation
{
  NONE,
  AVERAGE,
  HARMONIC
};

/// An input of a property given by variables: a number or a linear finite volume variable
struct Input
{
  /// Whether the input is a variable
  bool variable = false;
  /// The value of a number
  Real constant = 0;
  /// The cell values of a variable
  LinearFVCellValues values;
};

/// Where the face value of a variable input comes from on a range of faces
struct FaceSource
{
  enum Kind : int
  {
    /// The variable lives on neither side; as the host, which errors there, no value
    NONE,
    /// The linear interpolation of the two cell values
    INTERPOLATION,
    /// The cell value on the elem or neighbor side
    ELEM,
    NEIGHBOR,
    /// The value relation of a boundary condition, with the cell on the elem or neighbor side
    BC_ELEM,
    BC_NEIGHBOR
  };
  int kind = NONE;
  /// Index of the boundary condition in ExpressionData::bc_data
  int bc = -1;
  /// Index in the boundary condition's relation arrays of the first face of the range
  MOOSE_KOKKOS_INDEX_TYPE offset = 0;
};

/// The face lists of an object, for the face value sources
enum FaceList : unsigned int
{
  INTERNAL_FACES,
  BOUNDARY_FACES
};
}

namespace Moose::Kokkos
{

struct LinearFVPropertyInterface::ExpressionData
{
  /// The functor material with the expression, or nullptr for a single variable
  const KokkosParsedFunctorMaterial * material = nullptr;
  /// The inputs, in the order of the arguments of the expression
  Array<Input> inputs;
  /// The variables of the variable inputs, nullptr for numbers (host only)
  std::vector<const MooseLinearVariableFV<Real> *> variables;
  /// The Kokkos boundary conditions of the variable inputs with their compact index
  std::vector<LinearFVBoundaryCondition *> bcs;
  /// Value relation arrays of bcs
  Array<LinearFVBoundaryCondition::BoundaryData> bc_data;
  /// Face value sources of each (input, range) of the internal and boundary face lists
  Array2D<FaceSource> sources[2];
  /// The face lists the sources were built for
  const FaceRangeList * source_lists[2] = {nullptr, nullptr};
  /// Whether any input is a variable
  bool has_variables = false;
};

namespace
{
/// Evaluates a constant property
struct ConstantValue
{
  Real value;
  KOKKOS_FUNCTION Real cell(ContiguousElementID, const Real3 &) const { return value; }
  KOKKOS_FUNCTION Real
  face(ContiguousFaceID, unsigned int, unsigned int, MOOSE_KOKKOS_INDEX_TYPE, const Real3 &) const
  {
    return value;
  }
};

/// Evaluates a function property at a time
struct FunctionValue
{
  ReferenceWrapper<const KokkosParsedFunction> function;
  Real t;
  KOKKOS_FUNCTION Real cell(ContiguousElementID, const Real3 & point) const
  {
    return function->value(t, point);
  }
  KOKKOS_FUNCTION Real face(ContiguousFaceID,
                            unsigned int,
                            unsigned int,
                            MOOSE_KOKKOS_INDEX_TYPE,
                            const Real3 & point) const
  {
    return function->value(t, point);
  }
};

/// Evaluates a property given by variables, numbers and an expression at a time
struct ExpressionValue
{
  /// Whether there is an expression; a single input otherwise
  bool has_expression;
  RPNEvaluator evaluator;
  Real t;
  Array<System> systems;
  FVMesh mesh;
  Array<Input> inputs;
  Array2D<FaceSource> sources[2];
  Array<LinearFVBoundaryCondition::BoundaryData> bc_data;

  /// Get the cell value of an input, or NaN where the variable does not live
  KOKKOS_FUNCTION Real cellInput(const Input & input, const ContiguousElementID elem) const
  {
    if (!input.variable)
      return input.constant;
    const auto & system = systems[input.values.sys];
    const auto dof = system.getElemLocalDofIndex(elem, 0, input.values.var);
    if (dof == libMesh::DofObject::invalid_id)
      return std::numeric_limits<Real>::quiet_NaN();
    return system.getVectorDofValue(dof, input.values.tag);
  }

  /// Combine the input values
  KOKKOS_FUNCTION Real combine(const Real * const arguments, const Real3 & point) const
  {
    return has_expression ? evaluator.eval(t, point, 0, nullptr, arguments) : arguments[0];
  }

  KOKKOS_FUNCTION Real cell(const ContiguousElementID elem, const Real3 & point) const
  {
    Real arguments[KokkosParsedFunctorMaterial::max_functors];
    for (unsigned int i = 0; i < inputs.size(); ++i)
      arguments[i] = cellInput(inputs[i], elem);
    return combine(arguments, point);
  }

  KOKKOS_FUNCTION Real face(const ContiguousFaceID face,
                            const unsigned int list,
                            const unsigned int range,
                            const MOOSE_KOKKOS_INDEX_TYPE index,
                            const Real3 & point) const
  {
    Real arguments[KokkosParsedFunctorMaterial::max_functors];
    for (unsigned int i = 0; i < inputs.size(); ++i)
    {
      const auto & input = inputs[i];
      if (!input.variable)
      {
        arguments[i] = input.constant;
        continue;
      }

      // As MooseLinearVariableFV::evaluate() with a face argument
      const auto & source = sources[list](i, range);
      switch (source.kind)
      {
        case FaceSource::INTERPOLATION:
        {
          const Real gc = mesh.getGC(face);
          arguments[i] = gc * cellInput(input, mesh.getElem(face)) +
                         (1.0 - gc) * cellInput(input, mesh.getNeighbor(face));
          break;
        }
        case FaceSource::ELEM:
          arguments[i] = cellInput(input, mesh.getElem(face));
          break;
        case FaceSource::NEIGHBOR:
          arguments[i] = cellInput(input, mesh.getNeighbor(face));
          break;
        case FaceSource::BC_ELEM:
        case FaceSource::BC_NEIGHBOR:
        {
          const auto cell = source.kind == FaceSource::BC_ELEM ? mesh.getElem(face)
                                                               : mesh.getNeighbor(face);
          const auto & relation = bc_data[source.bc].value;
          const auto entry = source.offset + index;
          arguments[i] =
              relation.coefficient[entry] * cellInput(input, cell) + relation.source[entry];
          break;
        }
        default:
          arguments[i] = std::numeric_limits<Real>::quiet_NaN();
      }
    }
    return combine(arguments, point);
  }
};

/**
 * Fill the values of a property, without waiting for the work to finish
 * @param value The device evaluator of the property
 * @param fv_mesh The finite volume mesh data
 * @param cells The cell values to fill, or nullptr
 * @param faces The face values to fill, or nullptr
 * @param internal_faces The faces with the variable on both sides to fill, or nullptr
 * @param interpolation The interpolation of the cell values to the internal faces
 * @param boundary_faces The faces with the variable on one side to fill, or nullptr
 * @param offset_sign How each range of boundary_faces is evaluated: with the evaluation point
 * moved into the elem cell (-1), into the neighbor cell (1), at the face centroid (0), or as the
 * linear interpolation of the cell values (2)
 */
template <typename Value>
void
fill(const Value & value,
     const FVMesh & fv_mesh,
     Real * const cells,
     Real * const faces,
     const FaceRangeList * const internal_faces,
     const Interpolation interpolation,
     const FaceRangeList * const boundary_faces,
     const Array<int> & offset_sign)
{
  const auto mesh = fv_mesh;

  if (cells)
    ::Kokkos::parallel_for(
        ::Kokkos::RangePolicy<ExecSpace>(0, fv_mesh.getNumElements()),
        KOKKOS_LAMBDA(const ContiguousElementID elem) {
          cells[elem] = value.cell(elem, mesh.getElemCentroid(elem));
        });

  if (internal_faces && internal_faces->size())
  {
    const auto list = *internal_faces;
    ::Kokkos::parallel_for(
        ::Kokkos::RangePolicy<ExecSpace>(0, list.size()),
        KOKKOS_LAMBDA(const MOOSE_KOKKOS_INDEX_TYPE tid) {
          const auto range = list.range(tid);
          const auto face = list.face(tid, range);
          if (interpolation == Interpolation::NONE)
          {
            faces[face] = value.face(
                face, INTERNAL_FACES, range, list.indexInRange(tid, range), mesh.getCentroid(face));
            return;
          }

          // As FVGeometricAverage::interpolate() and FVHarmonicAverage::interpolate()
          const Real gc = mesh.getGC(face);
          const Real elem_value = cells[mesh.getElem(face)];
          const Real neighbor_value = cells[mesh.getNeighbor(face)];
          if (interpolation == Interpolation::AVERAGE)
          {
            faces[face] = gc * elem_value + (1.0 - gc) * neighbor_value;
            return;
          }

          const Real one_minus_gc = 1.0 - gc;
          // The smallest positive normalized value avoids dividing by zero and keeps the sign, as
          // on the host
          const auto safe = [](const Real v)
          {
            constexpr Real eps = std::numeric_limits<Real>::min();
            const Real magnitude = v < 0 ? -v : v;
            return ::Kokkos::copysign(magnitude > eps ? magnitude : eps, v);
          };
          const Real denom = gc / safe(elem_value) + one_minus_gc / safe(neighbor_value);
          faces[face] = 1.0 / denom;
        });
  }

  if (boundary_faces && boundary_faces->size())
  {
    const auto list = *boundary_faces;
    ::Kokkos::parallel_for(
        ::Kokkos::RangePolicy<ExecSpace>(0, list.size()),
        KOKKOS_LAMBDA(const MOOSE_KOKKOS_INDEX_TYPE tid) {
          const auto range = list.range(tid);
          const auto face = list.face(tid, range);
          if (offset_sign[range] == 2)
          {
            const Real gc = mesh.getGC(face);
            faces[face] =
                gc * cells[mesh.getElem(face)] + (1.0 - gc) * cells[mesh.getNeighbor(face)];
            return;
          }
          auto point = mesh.getCentroid(face);
          if (offset_sign[range])
          {
            // The offset of Function::evaluate() for a one-sided face argument between blocks:
            // 1e-8 along the normal, which points out of the elem cell
            auto offset = 1e-8 * mesh.getNormal(face);
            if (offset_sign[range] < 0)
              offset *= -1;
            point += offset;
          }
          faces[face] = value.face(face, BOUNDARY_FACES, range, list.indexInRange(tid, range), point);
        });
  }
}

/// Find the face range of the finite volume mesh starting at a face
unsigned int
meshRange(const FVMesh & fv_mesh, const ContiguousFaceID begin)
{
  const auto & ranges = fv_mesh.getFaceRanges();
  const auto it = std::upper_bound(ranges.begin(),
                                   ranges.end(),
                                   begin,
                                   [](const ContiguousFaceID f, const FaceRange & range)
                                   { return f < range.begin; });
  return std::distance(ranges.begin(), it) - 1;
}

/// Properties being updated, to detect a property that needs itself through a boundary condition
std::vector<const LinearFVPropertyInterface *> updating_properties;
}

LinearFVPropertyInterface::LinearFVPropertyInterface(const MooseObject & moose_object,
                                                     FEProblemBase & problem,
                                                     LinearFVVariableInterface & variables)
  : _lfvp_moose_object(moose_object), _lfvp_problem(problem), _lfvp_variables(&variables)
{
}

LinearFVPropertyInterface::LinearFVPropertyInterface(const LinearFVPropertyInterface & object)
  : _lfvp_moose_object(object._lfvp_moose_object),
    _lfvp_problem(object._lfvp_problem),
    _lfvp_variables(nullptr)
{
}

LinearFVPropertyInterface::~LinearFVPropertyInterface() = default;

void
LinearFVPropertyInterface::requestLinearFVProperty(const std::string & param,
                                                   LinearFVProperty & property,
                                                   const PropertyUse & use)
{
  if (use.interpolation && !dynamic_cast<const FVGeometricAverage *>(use.interpolation) &&
      !dynamic_cast<const FVHarmonicAverage *>(use.interpolation))
    _lfvp_moose_object.mooseError(
        "Kokkos linear finite volume properties are interpolated with FVGeometricAverage or "
        "FVHarmonicAverage only.");

  PropertyData data;
  data.param = param;
  data.handle = &property;
  data.use = use;

  const auto & name = _lfvp_moose_object.getParam<MooseFunctorName>(param);
  if (!MooseUtils::parsesToReal(name) && !_lfvp_problem.hasKokkosFunction(name))
  {
    // A variable or a functor material property, which may not exist yet
    data.deferred = true;
    _properties.push_back(std::move(data));
    return;
  }

  // Numbers become KokkosConstantFunction objects
  const auto & function = _lfvp_problem.getKokkosFunction<Moose::Kokkos::FunctionBase>(name);
  if (const auto * const constant = dynamic_cast<const KokkosConstantFunction *>(&function))
    property.constant = constant->value(0, Real3());
  else if (const auto * const parsed = dynamic_cast<const KokkosParsedFunction *>(&function))
    data.function = parsed;
  else
    _lfvp_moose_object.paramError(param,
                                  "The Kokkos function '",
                                  name,
                                  "' must be a KokkosParsedFunction or a KokkosConstantFunction.");

  // A constant needs no values unless it is interpolated, which the host does in floating point
  if (data.function || (use.interpolation && use.internal_faces))
    _properties.push_back(std::move(data));
}

void
LinearFVPropertyInterface::setupLinearFVProperties()
{
  for (auto & data : _properties)
    if (data.deferred)
    {
      setupExpression(data, _lfvp_moose_object.getParam<MooseFunctorName>(data.param));
      data.deferred = false;
    }
}

void
LinearFVPropertyInterface::setupExpression(PropertyData & data, const std::string & name)
{
  auto expression = std::make_shared<ExpressionData>();

  // The functors combined by the property
  std::vector<std::string> functors;
  if (_lfvp_problem.hasVariable(name))
    functors.push_back(name);
  else
  {
    for (const auto & material : _lfvp_problem.getMaterialWarehouse().getObjects())
    {
      const auto * const parsed = dynamic_cast<const KokkosParsedFunctorMaterial *>(material.get());
      // Functor materials are also built for the displaced problem, which defines the same names
      if (!parsed || parsed->propertyName() != name ||
          parsed->getParam<SubProblem *>("_subproblem") != &_lfvp_problem)
        continue;
      if (expression->material)
        _lfvp_moose_object.paramError(data.param,
                                      "The property '",
                                      name,
                                      "' is defined by several KokkosParsedFunctorMaterial objects, "
                                      "which Kokkos linear finite volume objects do not support.");
      expression->material = parsed;
    }
    if (!expression->material)
      _lfvp_moose_object.paramError(
          data.param,
          "'",
          name,
          "' is not a number, a Kokkos function, a linear finite volume variable or the "
          "property_name of a KokkosParsedFunctorMaterial.");
    functors = expression->material->functorNames();
  }

  // The inputs and the Kokkos boundary conditions of the variable inputs
  std::vector<Input> inputs(functors.size());
  expression->variables.resize(functors.size(), nullptr);
  for (const auto i : index_range(functors))
  {
    Real constant;
    if (MooseUtils::parsesToReal(functors[i], &constant))
    {
      inputs[i].constant = constant;
      continue;
    }

    const auto * const variable =
        _lfvp_problem.hasVariable(functors[i])
            ? dynamic_cast<const MooseLinearVariableFV<Real> *>(
                  &_lfvp_problem.getVariable(0, functors[i]))
            : nullptr;
    if (!variable)
      _lfvp_moose_object.paramError(
          data.param,
          "The functor '",
          functors[i],
          "' of the property '",
          name,
          "' must be a number or a linear finite volume variable (MooseLinearVariableFVReal).");

    inputs[i].variable = true;
    inputs[i].values = _lfvp_variables->getLinearFVCellValues(functors[i]);
    expression->variables[i] = variable;
    expression->has_variables = true;

    std::vector<::LinearFVBoundaryCondition *> host_bcs;
    _lfvp_problem.theWarehouse()
        .query()
        .template condition<AttribSysNum>(variable->sys().number())
        .template condition<AttribVar>(variable->number())
        .template condition<AttribSystem>("LinearFVBoundaryCondition")
        .template condition<AttribKokkos>(false)
        .template condition<AttribThread>(0)
        .queryInto(host_bcs);
    for (const auto * const bc : host_bcs)
      _lfvp_moose_object.paramError(
            data.param,
            "The variable '",
            functors[i],
            "' of the property '",
            name,
            "' has the host boundary condition '",
            bc->name(),
            "', whose boundary values Kokkos linear finite volume objects cannot read.");
    // The Kokkos boundary conditions are found by system and filtered by variable, as for the
    // Green-Gauss gradients
    std::vector<LinearFVBoundaryCondition *> kokkos_bcs;
    _lfvp_problem.theWarehouse()
        .query()
        .template condition<AttribSysNum>(variable->sys().number())
        .template condition<AttribSystem>("LinearFVBoundaryCondition")
        .template condition<AttribKokkos>(true)
        .template condition<AttribThread>(0)
        .queryInto(kokkos_bcs);
    for (auto * const bc : kokkos_bcs)
    {
      if (bc->variable().sys() != variable->sys().number() ||
          bc->variable().var() != variable->number())
        continue;
      if (!bc->hasBoundaryValue())
        _lfvp_moose_object.paramError(data.param,
                                      "The boundary condition '",
                                      bc->name(),
                                      "' of the variable '",
                                      functors[i],
                                      "' does not provide boundary values, which the property '",
                                      name,
                                      "' needs.");
      expression->bcs.push_back(bc);
    }
  }
  expression->inputs = inputs;

  data.expression = expression;
}

void
LinearFVPropertyInterface::updateLinearFVProperties(const Real t,
                                                    const FaceRangeList * const internal_faces,
                                                    const FaceRangeList * const boundary_faces,
                                                    const std::vector<bool> & boundary_flips)
{
  if (_properties.empty())
    return;

  const auto & fv_mesh = _lfvp_problem.getKokkosFVMesh();
  const auto num_elems = fv_mesh.getNumElements();
  const auto num_faces = fv_mesh.getNumFaces();

  // The boundary conditions whose boundary values the properties read, updated before the fill
  std::vector<LinearFVBoundaryCondition *> updated_bcs;

  for (auto & data : _properties)
  {
    mooseAssert(!data.deferred, "setupLinearFVProperties() must be called first");

    // The host evaluates a functor material on a face with the variable on both sides, and on a
    // face between blocks without a side, as the linear interpolation of its cell values
    // (PiecewiseByBlockLambdaFunctor::evaluate())
    const bool interpolates_cells = data.expression && data.expression->material;
    const bool needs_cells =
        data.use.cells ||
        ((data.use.interpolation || interpolates_cells) && data.use.internal_faces) ||
        (interpolates_cells && data.use.boundary_faces && data.use.unsided_boundary_faces);
    const bool needs_faces = data.use.internal_faces || data.use.boundary_faces;
    // Constants are filled once, functions when the time changes, variables every time, and all
    // when the mesh changes
    const bool resized = (needs_cells && data.cell_values.size() != num_elems) ||
                         (needs_faces && data.face_values.size() != num_faces);
    const bool time_dependent = data.function || data.expression;
    const bool solution_dependent = data.expression && data.expression->has_variables;
    if (data.filled && !resized && !solution_dependent && (!time_dependent || data.time == t))
      continue;

    if (needs_cells && data.cell_values.size() != num_elems)
      data.cell_values.createDevice(num_elems);
    if (needs_faces && data.face_values.size() != num_faces)
      data.face_values.createDevice(num_faces);
    data.handle->cell_values = needs_cells ? data.cell_values.deviceData() : nullptr;
    data.handle->face_values =
        data.face_values.isAlloc() ? data.face_values.deviceData() : nullptr;

    Real * const cells = needs_cells ? data.cell_values.deviceData() : nullptr;
    Real * const faces = data.face_values.isAlloc() ? data.face_values.deviceData() : nullptr;
    const auto * const internal = data.use.internal_faces ? internal_faces : nullptr;
    const auto * const boundary = data.use.boundary_faces ? boundary_faces : nullptr;

    auto interpolation = Interpolation::NONE;
    if (dynamic_cast<const FVGeometricAverage *>(data.use.interpolation) || (!data.use.interpolation && interpolates_cells))
      interpolation = Interpolation::AVERAGE;
    else if (dynamic_cast<const FVHarmonicAverage *>(data.use.interpolation))
      interpolation = Interpolation::HARMONIC;

    // Direction of the offset of the evaluation point of a function on each range of boundary
    // faces: on faces between blocks, into the cell of the variable, or none
    auto & offset_sign = data.offset_sign;
    if (boundary && boundary->numRanges())
    {
      const auto & ranges = fv_mesh.getFaceRanges();
      if (offset_sign.size() != boundary->numRanges())
        offset_sign.create(boundary->numRanges());
      for (const auto r : make_range(boundary->numRanges()))
      {
        const auto & range = ranges[meshRange(fv_mesh, boundary->begin(r))];
        const bool between_blocks = range.neighbor_subdomain != Moose::INVALID_BLOCK_ID &&
                                    range.neighbor_subdomain != range.elem_subdomain;
        if (between_blocks && interpolates_cells && data.use.unsided_boundary_faces)
          offset_sign[r] = 2;
        else
          offset_sign[r] = between_blocks && data.function && !data.use.unsided_boundary_faces
                               ? (boundary_flips[r] ? 1 : -1)
                               : 0;
      }
      offset_sign.copyToDevice();
    }

    if (data.function)
      fill(FunctionValue{*data.function, t},
           fv_mesh,
           cells,
           faces,
           internal,
           interpolation,
           boundary,
           offset_sign);
    else if (data.expression)
    {
      auto & expression = *data.expression;

      // The face value sources of the variable inputs on the ranges of each face list
      const FaceRangeList * const lists[2] = {internal, boundary};
      for (const auto l : make_range(2))
      {
        const auto * const list = lists[l];
        if (!list || !list->numRanges() || expression.source_lists[l] == list)
          continue;
        expression.source_lists[l] = list;

        auto & sources = expression.sources[l];
        sources.create(expression.inputs.size(), list->numRanges());
        for (const auto i : index_range(expression.variables))
        {
          const auto * const variable = expression.variables[i];
          for (const auto r : make_range(list->numRanges()))
          {
            sources(i, r) = FaceSource();
            if (!variable)
              continue;

            const auto mesh_range = meshRange(fv_mesh, list->begin(r));
            const auto & range = fv_mesh.getFaceRanges()[mesh_range];
            const bool on_elem = variable->activeOnSubdomain(range.elem_subdomain);
            const bool on_neighbor = range.neighbor_subdomain != Moose::INVALID_BLOCK_ID &&
                                     variable->activeOnSubdomain(range.neighbor_subdomain);
            FaceSource source;
            if (on_elem && on_neighbor)
              source.kind = FaceSource::INTERPOLATION;
            else if (on_elem || on_neighbor)
            {
              source.kind = on_elem ? FaceSource::ELEM : FaceSource::NEIGHBOR;
              for (const auto b : index_range(expression.bcs))
              {
                auto * const bc = expression.bcs[b];
                if (bc->variable().sys() != variable->sys().number() ||
                    bc->variable().var() != variable->number())
                  continue;
                bc->setupFaces();
                if (const auto offset = bc->faceRangeOffset(mesh_range))
                {
                  source.kind = on_elem ? FaceSource::BC_ELEM : FaceSource::BC_NEIGHBOR;
                  source.bc = b;
                  source.offset = *offset;
                }
              }
            }
            sources(i, r) = source;
          }
        }
        sources.copyToDevice();
      }

      // The boundary values the variable inputs read, with the current solution
      if (expression.bc_data.size() != expression.bcs.size())
        expression.bc_data.create(expression.bcs.size());
      for (const auto b : index_range(expression.bcs))
      {
        auto * const bc = expression.bcs[b];
        if (std::find(updating_properties.begin(), updating_properties.end(), bc) !=
            updating_properties.end())
          _lfvp_moose_object.paramError(data.param,
                                        "The property reads the boundary values of the boundary "
                                        "condition '",
                                        bc->name(),
                                        "', which need the property itself.");
        if (std::find(updated_bcs.begin(), updated_bcs.end(), bc) == updated_bcs.end())
        {
          updating_properties.push_back(this);
          bc->acquireGradients();
          bc->computeBoundaryValueData();
          updating_properties.pop_back();
          updated_bcs.push_back(bc);
        }
        expression.bc_data[b] = bc->boundaryData();
      }
      expression.bc_data.copyToDevice();

      ExpressionValue value{static_cast<bool>(expression.material),
                            {},
                            t,
                            _lfvp_problem.getKokkosSystems(),
                            fv_mesh,
                            expression.inputs,
                            {expression.sources[0], expression.sources[1]},
                            expression.bc_data};
      if (expression.material)
        value.evaluator = expression.material->evaluator();
      fill(value, fv_mesh, cells, faces, internal, interpolation, boundary, offset_sign);
    }
    else
      fill(ConstantValue{data.handle->constant},
           fv_mesh,
           cells,
           faces,
           internal,
           interpolation,
           boundary,
           offset_sign);

    data.time = t;
    data.filled = true;
  }

  // The boundary conditions read their gradients until the fill finished
  if (!updated_bcs.empty())
  {
    ::Kokkos::fence();
    for (auto * const bc : updated_bcs)
      bc->releaseGradients();
  }
}

} // namespace Moose::Kokkos
