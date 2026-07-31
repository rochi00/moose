//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LevelSetInterfaceMarker.h"

#include "FVUtils.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "MassConservedLevelSetCorrector.h"
#include "SubProblem.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_object.h"
#include "libmesh/numeric_vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

registerMooseObject("NavierStokesApp", LevelSetInterfaceMarker);

InputParameters
LevelSetInterfaceMarker::validParams()
{
  InputParameters params = Marker::validParams();
  params.addClassDescription(
      "Refines a local-cell-width band around a linear-FV signed-distance interface and coarsens "
      "outside a wider hysteresis band.");
  params.addRequiredParam<VariableName>("level_set_variable", "Linear-FV signed-distance field.");
  params.addRangeCheckedParam<Real>("refine_band_cells",
                                    4.0,
                                    "refine_band_cells>0",
                                    "Half-width of the refinement band in local cell widths.");
  params.addRangeCheckedParam<Real>(
      "coarsen_band_cells",
      7.0,
      "coarsen_band_cells>0",
      "Distance beyond which cells may coarsen, in local cell widths.");
  params.addParam<bool>(
      "conservative_prolongation",
      true,
      "Use a bounded linear reconstruction that preserves the parent-cell average when the "
      "cell-centered level set is prolonged to newly refined Cartesian children.");
  params.addRangeCheckedParam<Real>(
      "conservative_interface_width",
      0.0,
      "conservative_interface_width>=0",
      "Half-width of the regularized Heaviside whose quadrature cell average is conserved during "
      "transient refinement and coarsening. A value of zero retains signed-distance-only affine "
      "prolongation.");
  params.addParam<UserObjectName>(
      "level_set_corrector",
      "",
      "Mass-conserved signed-distance corrector applied to every BDF2 state immediately after a "
      "conservative material-fraction mesh transfer.");
  params.addRelationshipManager("ElementSideNeighborLayers",
                                Moose::RelationshipManagerType::GEOMETRIC |
                                    Moose::RelationshipManagerType::ALGEBRAIC |
                                    Moose::RelationshipManagerType::COUPLING,
                                [](const InputParameters &, InputParameters & rm_params)
                                { rm_params.set<unsigned short>("layers") = 1; });
  return params;
}

LevelSetInterfaceMarker::LevelSetInterfaceMarker(const InputParameters & parameters)
  : Marker(parameters),
    _level_set_variable(dynamic_cast<MooseLinearVariableFVReal &>(
        _subproblem.getVariable(_tid, getParam<VariableName>("level_set_variable")))),
    _level_set_system(dynamic_cast<LinearSystem &>(_level_set_variable.sys())),
    _system_number(_level_set_variable.sys().number()),
    _variable_number(_level_set_variable.number()),
    _refine_band_cells(getParam<Real>("refine_band_cells")),
    _coarsen_band_cells(getParam<Real>("coarsen_band_cells")),
    _conservative_prolongation(getParam<bool>("conservative_prolongation")),
    _conservative_interface_width(getParam<Real>("conservative_interface_width")),
    _level_set_corrector(
        !getParam<UserObjectName>("level_set_corrector").empty()
            ? const_cast<MassConservedLevelSetCorrector *>(
                  &getUserObject<MassConservedLevelSetCorrector>("level_set_corrector"))
            : nullptr)
{
  if (_coarsen_band_cells <= _refine_band_cells)
    paramError("coarsen_band_cells",
               "The coarsening band must be wider than the refinement band to provide "
               "refinement/coarsening hysteresis.");
  if (_conservative_interface_width > 0.0 && !_level_set_corrector)
    paramError("level_set_corrector",
               "A corrector is required when conservative material-fraction transfer is enabled "
               "so every BDF2 state is reconstructed as signed distance before the next solve.");

  addMooseVariableDependency(&_level_set_variable);
}

void
LevelSetInterfaceMarker::initialSetup()
{
  initializeProlongationVectors();
}

void
LevelSetInterfaceMarker::initializeProlongationVectors()
{
  if (!_conservative_prolongation)
    return;

  const unsigned int dimension = _mesh.dimension();
  for (unsigned int state = _prolongation_gradients.size();
       _level_set_system.hasSolutionState(state);
       ++state)
  {
    _prolongation_gradients.emplace_back(dimension);
    for (const auto axis : make_range(dimension))
      _prolongation_gradients[state][axis] =
          &_level_set_system.addVector("level_set_amr_gradient_" + name() + "_state_" +
                                           std::to_string(state) + "_axis_" + std::to_string(axis),
                                       true,
                                       GHOSTED);
    _projected_material_fractions.push_back(&_level_set_system.addVector(
        "level_set_amr_fraction_" + name() + "_state_" + std::to_string(state), true, GHOSTED));
  }
}

void
LevelSetInterfaceMarker::markerSetup()
{
  if (_tid == 0 && _conservative_prolongation)
  {
    initializeProlongationVectors();
    cacheProlongationGradients();
  }
}

void
LevelSetInterfaceMarker::cacheProlongationGradients()
{
  const unsigned int dimension = _mesh.dimension();
  if (dimension < 1 || dimension > 3)
    mooseError(name(), ": conservative FV prolongation requires dimension one, two, or three.");

  for (const auto & state_vectors : _prolongation_gradients)
    for (auto * const vector : state_vectors)
      vector->zero();
  for (auto * const vector : _projected_material_fractions)
    vector->zero();

  for (const auto & elem_info : _mesh.elemInfoVector())
  {
    const auto * const elem = elem_info->elem();
    if (!_level_set_variable.hasBlocks(elem->subdomain_id()) ||
        elem->processor_id() != processor_id())
      continue;
    if (_subproblem.getCoordSystem(elem->subdomain_id()) != Moose::COORD_XYZ)
      mooseError(name(),
                 ": conservative level-set prolongation currently requires Cartesian XYZ "
                 "coordinates.");

    const dof_id_type dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (dof == DofObject::invalid_id)
      continue;

    Real scale = elem->hmin();
    std::array<Real, LIBMESH_DIM> cell_extent = {};
    for (const auto axis : make_range(dimension))
    {
      Real minimum = std::numeric_limits<Real>::max();
      Real maximum = std::numeric_limits<Real>::lowest();
      for (const auto node : make_range(elem->n_nodes()))
      {
        minimum = std::min(minimum, elem->point(node)(axis));
        maximum = std::max(maximum, elem->point(node)(axis));
      }
      cell_extent[axis] = maximum - minimum;
      if (cell_extent[axis] <= 0.0)
        mooseError(name(), ": zero Cartesian extent on element ", elem->id(), ".");
    }

    for (const auto state : index_range(_prolongation_gradients))
    {
      const auto state_arg =
          Moose::StateArg{cast_int<unsigned int>(state), Moose::SolutionIterationType::Time};
      const Real center_value = _level_set_variable.getElemValue(*elem_info, state_arg);
      Real minimum_value = center_value;
      Real maximum_value = center_value;
      DenseMatrix<Real> normal_matrix(dimension, dimension);
      DenseVector<Real> right_hand_side(dimension);
      std::array<bool, 2 * LIBMESH_DIM> has_neighbor = {};

      auto inspect_face = [this,
                           dimension,
                           scale,
                           elem_info,
                           state_arg,
                           center_value,
                           &normal_matrix,
                           &right_hand_side,
                           &has_neighbor,
                           &minimum_value,
                           &maximum_value](const Elem &,
                                           const Elem * const neighbor,
                                           const FaceInfo *,
                                           const Point & surface_vector,
                                           const Real,
                                           const bool)
      {
        unsigned int axis = 0;
        for (const auto component : make_range(1u, dimension))
          if (std::abs(surface_vector(component)) > std::abs(surface_vector(axis)))
            axis = component;
        for (const auto component : make_range(dimension))
          if (component != axis &&
              std::abs(surface_vector(component)) > 1e-10 * std::abs(surface_vector(axis)))
            mooseError(name(),
                       ": conservative level-set prolongation requires axis-aligned Cartesian "
                       "faces, but element ",
                       elem_info->elem()->id(),
                       " has surface vector ",
                       surface_vector,
                       ".");

        const bool positive = surface_vector(axis) > 0.0;
        if (!neighbor || !_level_set_variable.hasBlocks(neighbor->subdomain_id()))
          return;
        has_neighbor[2 * axis + (positive ? 1 : 0)] = true;

        const auto & neighbor_info = _mesh.elemInfo(neighbor->id());
        const Real neighbor_value = _level_set_variable.getElemValue(neighbor_info, state_arg);
        minimum_value = std::min(minimum_value, neighbor_value);
        maximum_value = std::max(maximum_value, neighbor_value);

        const Point normalized_offset = (neighbor_info.centroid() - elem_info->centroid()) / scale;
        const Real weight = surface_vector.norm() / std::max(normalized_offset.norm_sq(),
                                                             std::numeric_limits<Real>::epsilon());
        for (const auto row : make_range(dimension))
        {
          right_hand_side(row) += weight * normalized_offset(row) * (neighbor_value - center_value);
          for (const auto column : make_range(dimension))
            normal_matrix(row, column) +=
                weight * normalized_offset(row) * normalized_offset(column);
        }
      };

      Moose::FV::loopOverElemFaceInfo(
          *elem, _mesh, inspect_face, Moose::COORD_XYZ, libMesh::invalid_uint);

      DenseVector<Real> scaled_gradient(dimension);
      normal_matrix.svd_solve(right_hand_side, scaled_gradient);
      RealVectorValue gradient;
      for (const auto axis : make_range(dimension))
        gradient(axis) = scaled_gradient(axis) / scale;

      // Supply linearly extrapolated ghost-cell bounds at physical boundaries. This retains exact
      // linear reconstruction there instead of incorrectly treating a boundary cell as an
      // extremum merely because it has neighbors on only one side.
      for (const auto axis : make_range(dimension))
        for (const bool positive : {false, true})
          if (!has_neighbor[2 * axis + (positive ? 1 : 0)])
          {
            const Real ghost_value =
                center_value + (positive ? 1.0 : -1.0) * gradient(axis) * cell_extent[axis];
            minimum_value = std::min(minimum_value, ghost_value);
            maximum_value = std::max(maximum_value, ghost_value);
          }

      Real maximum_child_delta = 0.0;
      for (const auto axis : make_range(dimension))
        maximum_child_delta += 0.25 * std::abs(gradient(axis)) * cell_extent[axis];
      Real limiter = 1.0;
      if (maximum_child_delta > std::numeric_limits<Real>::epsilon())
        limiter = std::min({limiter,
                            std::max(0.0, maximum_value - center_value) / maximum_child_delta,
                            std::max(0.0, center_value - minimum_value) / maximum_child_delta});

      for (const auto axis : make_range(dimension))
        _prolongation_gradients[state][axis]->set(dof, limiter * gradient(axis));
      _projected_material_fractions[state]->set(
          dof,
          cellAverageMaterialFraction(center_value, limiter * gradient, cell_extent, dimension));
    }
  }

  for (const auto & state_vectors : _prolongation_gradients)
    for (auto * const vector : state_vectors)
      vector->close();
  for (auto * const vector : _projected_material_fractions)
    vector->close();
}

void
LevelSetInterfaceMarker::meshChanged()
{
  if (_tid != 0 || !_conservative_prolongation)
    return;

  prolongRefinedCells();
  restrictCoarsenedCells();
  if (_conservative_interface_width > 0.0 && _fe_problem.timeStep() > 0)
    _level_set_corrector->correctAfterMeshChange();
  cacheProlongationGradients();
}

Real
LevelSetInterfaceMarker::regularizedHeaviside(const Real phi) const
{
  if (_conservative_interface_width == 0.0)
    return phi >= 0.0 ? 1.0 : 0.0;
  if (phi <= -_conservative_interface_width)
    return 0.0;
  if (phi >= _conservative_interface_width)
    return 1.0;

  const Real normalized_phi = phi / _conservative_interface_width;
  return 0.5 * (1.0 + normalized_phi + std::sin(libMesh::pi * normalized_phi) / libMesh::pi);
}

Real
LevelSetInterfaceMarker::conservativeRefinementShift(const std::vector<Real> & child_values,
                                                     const Real parent_fraction) const
{
  mooseAssert(!child_values.empty(), "At least one child value is required for AMR transfer.");
  const auto [minimum, maximum] = std::minmax_element(child_values.begin(), child_values.end());
  // A saturated indicator contains no subcell interface location from which a signed distance can
  // be reconstructed. Retain the affine level-set reconstruction in pure cells; only mixed cells
  // possess enough material information to define a conservative local correction.
  if (parent_fraction <= 0.0 || parent_fraction >= 1.0)
    return 0.0;

  auto mean_fraction = [this, &child_values](const Real shift)
  {
    Real sum = 0.0;
    for (const Real value : child_values)
      sum += regularizedHeaviside(value + shift);
    return sum / child_values.size();
  };

  Real lower = -_conservative_interface_width - *maximum;
  Real upper = _conservative_interface_width - *minimum;
  for (const auto iteration : make_range(80))
  {
    const Real midpoint = 0.5 * (lower + upper);
    if (mean_fraction(midpoint) < parent_fraction)
      lower = midpoint;
    else
      upper = midpoint;
    (void)iteration;
  }
  return 0.5 * (lower + upper);
}

Real
LevelSetInterfaceMarker::cellAverageMaterialFraction(
    const Real center_value,
    const RealVectorValue & gradient,
    const std::array<Real, LIBMESH_DIM> & cell_extent,
    const unsigned int dimension) const
{
  Real fraction_sum = 0.0;
  const unsigned int number_of_children = 1u << dimension;
  const unsigned int points_per_child = 1u << dimension;
  for (const auto child_index : make_range(number_of_children))
  {
    Point child_center_offset;
    for (const auto axis : make_range(dimension))
      child_center_offset(axis) = (child_index & (1u << axis) ? 0.25 : -0.25) * cell_extent[axis];

    for (const auto quadrature_index : make_range(points_per_child))
    {
      Real quadrature_value = center_value;
      for (const auto axis : make_range(dimension))
      {
        const Real quadrature_offset = (quadrature_index & (1u << axis) ? 1.0 : -1.0) *
                                       cell_extent[axis] / (4.0 * std::sqrt(3.0));
        quadrature_value += gradient(axis) * (child_center_offset(axis) + quadrature_offset);
      }
      fraction_sum += regularizedHeaviside(quadrature_value);
    }
  }
  return fraction_sum / (number_of_children * points_per_child);
}

void
LevelSetInterfaceMarker::prolongRefinedCells()
{
  const auto * const refined_elements = _mesh.refinedElementRange();
  if (!refined_elements)
    return;

  const unsigned int dimension = _mesh.dimension();
  for (const auto state : index_range(_prolongation_gradients))
  {
    auto & solution = _level_set_system.solutionState(state);
    for (const auto * const parent : *refined_elements)
    {
      if (!_level_set_variable.hasBlocks(parent->subdomain_id()))
        continue;
      const Point parent_centroid = parent->vertex_average();
      Real conservative_shift = 0.0;
      if (_conservative_interface_width > 0.0 && _fe_problem.timeStep() > 0)
      {
        const Elem * representative_child = nullptr;
        for (const auto child_index : make_range(parent->n_children()))
        {
          const auto * const child = parent->child_ptr(child_index);
          if (child && child != remote_elem && child->active() &&
              child->processor_id() == processor_id())
          {
            representative_child = child;
            break;
          }
        }
        if (!representative_child)
          continue;

        const auto & representative_info = _mesh.elemInfo(representative_child->id());
        const dof_id_type representative_dof =
            representative_info.dofIndices()[_system_number][_variable_number];
        const Real parent_value = solution(representative_dof);
        const Real parent_fraction = (*_projected_material_fractions[state])(representative_dof);

        std::array<Real, LIBMESH_DIM> parent_extent = {};
        for (const auto axis : make_range(dimension))
        {
          Real minimum = std::numeric_limits<Real>::max();
          Real maximum = std::numeric_limits<Real>::lowest();
          for (const auto node : make_range(parent->n_nodes()))
          {
            minimum = std::min(minimum, parent->point(node)(axis));
            maximum = std::max(maximum, parent->point(node)(axis));
          }
          parent_extent[axis] = maximum - minimum;
        }

        std::vector<Real> child_values;
        child_values.reserve(1u << (2 * dimension));
        for (const auto child_index : make_range(1u << dimension))
        {
          Point child_center_offset;
          for (const auto axis : make_range(dimension))
            child_center_offset(axis) =
                (child_index & (1u << axis) ? 0.25 : -0.25) * parent_extent[axis];

          for (const auto quadrature_index : make_range(1u << dimension))
          {
            Real child_value = parent_value;
            for (const auto axis : make_range(dimension))
            {
              const Real quadrature_offset = (quadrature_index & (1u << axis) ? 1.0 : -1.0) *
                                             parent_extent[axis] / (4.0 * std::sqrt(3.0));
              child_value += (*_prolongation_gradients[state][axis])(representative_dof) *
                             (child_center_offset(axis) + quadrature_offset);
            }
            child_values.push_back(child_value);
          }
        }
        conservative_shift = conservativeRefinementShift(child_values, parent_fraction);
      }

      for (const auto child_index : make_range(parent->n_children()))
      {
        const auto * const child = parent->child_ptr(child_index);
        if (!child || child == remote_elem || !child->active() ||
            child->processor_id() != processor_id())
          continue;

        const auto & child_info = _mesh.elemInfo(child->id());
        const dof_id_type dof = child_info.dofIndices()[_system_number][_variable_number];
        const Point offset = child_info.centroid() - parent_centroid;
        Real reconstructed_value = solution(dof);
        for (const auto axis : make_range(dimension))
          reconstructed_value += (*_prolongation_gradients[state][axis])(dof)*offset(axis);
        reconstructed_value += conservative_shift;
        solution.set(dof, reconstructed_value);
      }
    }
    solution.close();
  }

  auto & current_local_solution = *_level_set_system.system().current_local_solution;
  _level_set_system.solution().localize(current_local_solution,
                                        _level_set_system.dofMap().get_send_list());
  _level_set_system.setSolution(current_local_solution);
}

void
LevelSetInterfaceMarker::restrictCoarsenedCells()
{
  if (_conservative_interface_width == 0.0 || _fe_problem.timeStep() == 0)
    return;

  const auto * const coarsened_elements = _mesh.coarsenedElementRange();
  if (!coarsened_elements)
    return;

  for (const auto state : index_range(_projected_material_fractions))
  {
    auto & solution = _level_set_system.solutionState(state);
    for (const auto * const parent : *coarsened_elements)
    {
      if (!_level_set_variable.hasBlocks(parent->subdomain_id()) ||
          parent->processor_id() != processor_id())
        continue;

      const auto & parent_info = _mesh.elemInfo(parent->id());
      const dof_id_type dof = parent_info.dofIndices()[_system_number][_variable_number];
      const Real projected_value = solution(dof);
      const Real parent_fraction = (*_projected_material_fractions[state])(dof);
      const unsigned int dimension = _mesh.dimension();
      std::array<Real, LIBMESH_DIM> parent_extent = {};
      for (const auto axis : make_range(dimension))
      {
        Real minimum = std::numeric_limits<Real>::max();
        Real maximum = std::numeric_limits<Real>::lowest();
        for (const auto node : make_range(parent->n_nodes()))
        {
          minimum = std::min(minimum, parent->point(node)(axis));
          maximum = std::max(maximum, parent->point(node)(axis));
        }
        parent_extent[axis] = maximum - minimum;
      }

      std::vector<Real> quadrature_values;
      quadrature_values.reserve(1u << dimension);
      for (const auto quadrature_index : make_range(1u << dimension))
      {
        Real quadrature_value = projected_value;
        for (const auto axis : make_range(dimension))
        {
          const Real quadrature_offset = (quadrature_index & (1u << axis) ? 1.0 : -1.0) *
                                         parent_extent[axis] / (2.0 * std::sqrt(3.0));
          quadrature_value += (*_prolongation_gradients[state][axis])(dof)*quadrature_offset;
        }
        quadrature_values.push_back(quadrature_value);
      }
      solution.set(
          dof, projected_value + conservativeRefinementShift(quadrature_values, parent_fraction));
    }
    solution.close();
  }

  auto & current_local_solution = *_level_set_system.system().current_local_solution;
  _level_set_system.solution().localize(current_local_solution,
                                        _level_set_system.dofMap().get_send_list());
  _level_set_system.setSolution(current_local_solution);
}

Marker::MarkerValue
LevelSetInterfaceMarker::computeElementMarker()
{
  const auto state = Moose::currentState();
  const auto & elem_info = _mesh.elemInfo(_current_elem->id());
  const Real phi = _level_set_variable.getElemValue(elem_info, state);
  if (!std::isfinite(phi))
    mooseError(name(), ": non-finite level-set value in element ", _current_elem->id(), ".");

  bool crosses_interface = phi == 0.0;
  auto inspect_face = [this, state, phi, &crosses_interface](const Elem &,
                                                             const Elem * const neighbor,
                                                             const FaceInfo *,
                                                             const Point &,
                                                             const Real,
                                                             const bool)
  {
    if (!neighbor || !_level_set_variable.hasBlocks(neighbor->subdomain_id()))
      return;

    const auto & neighbor_info = _mesh.elemInfo(neighbor->id());
    const Real neighbor_phi = _level_set_variable.getElemValue(neighbor_info, state);
    if (!std::isfinite(neighbor_phi))
      mooseError(name(), ": non-finite level-set value in element ", neighbor->id(), ".");

    crosses_interface = crosses_interface || neighbor_phi == 0.0 || phi * neighbor_phi < 0.0;
  };

  const auto coord_type = _subproblem.getCoordSystem(_current_elem->subdomain_id());
  Moose::FV::loopOverElemFaceInfo(*_current_elem,
                                  _mesh,
                                  inspect_face,
                                  coord_type,
                                  coord_type == Moose::COORD_RZ
                                      ? _subproblem.getAxisymmetricRadialCoord()
                                      : libMesh::invalid_uint);

  if (crosses_interface)
    return REFINE;

  const Real local_cell_width = _current_elem->hmin();
  const Real interface_distance = std::abs(phi);
  if (interface_distance <= _refine_band_cells * local_cell_width)
    return REFINE;
  if (interface_distance >= _coarsen_band_cells * local_cell_width)
    return COARSEN;
  return DO_NOTHING;
}
