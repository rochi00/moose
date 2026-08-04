//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MassConservedLevelSetCorrector.h"

#include "LinearSystem.h"
#include "MooseLinearVariableFV.h"
#include "FVUtils.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/remote_elem.h"

registerMooseObject("NavierStokesApp", MassConservedLevelSetCorrector);

namespace
{
Real
minmod(const Real first, const Real second)
{
  if (first * second <= 0.0)
    return 0.0;
  return std::abs(first) < std::abs(second) ? first : second;
}

int
strictSign(const Real value)
{
  return (value > 0.0) - (value < 0.0);
}

unsigned int
numberOfQuadraticCoefficients(const unsigned int dimension)
{
  return dimension * (dimension + 3) / 2;
}

std::vector<Real>
quadraticBasis(const Point & offset, const Real scale, const unsigned int dimension)
{
  std::vector<Real> basis(numberOfQuadraticCoefficients(dimension));
  std::array<Real, LIBMESH_DIM> coordinate = {};
  for (const auto axis : make_range(dimension))
  {
    coordinate[axis] = offset(axis) / scale;
    basis[axis] = coordinate[axis];
    basis[dimension + axis] = 0.5 * coordinate[axis] * coordinate[axis];
  }

  unsigned int coefficient = 2 * dimension;
  for (const auto first : make_range(dimension))
    for (const auto second : make_range(first + 1, dimension))
      basis[coefficient++] = coordinate[first] * coordinate[second];
  return basis;
}

unsigned int
directionIndex(const unsigned int axis, const bool positive)
{
  return 2 * axis + (positive ? 1 : 0);
}
}

InputParameters
MassConservedLevelSetCorrector::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params += BlockRestrictable::validParams();
  params.addClassDescription(
      "Conserves the material mass represented by a signed-distance level set through a uniform "
      "interface shift.");
  params.addRequiredParam<VariableName>("level_set_variable",
                                        "Linear finite-volume signed-distance variable.");
  params.addRequiredParam<MooseFunctorName>(
      "material_fraction", "Name of the regularized-Heaviside material-fraction functor.");
  params.addRequiredParam<MooseFunctorName>(
      "material_density", "Solid/liquid material density used in the mass constraint.");
  params.addRequiredParam<MooseFunctorName>(
      "volumetric_face_flux", "Face flux used by the level-set material-advection equation.");
  params.addRangeCheckedParam<Real>(
      "interface_width", "interface_width>0", "Half-width of the regularized Heaviside.");
  params.addRangeCheckedParam<Real>("relative_tolerance",
                                    1e-12,
                                    "relative_tolerance>0",
                                    "Relative tolerance for the material-mass constraint.");
  params.addRangeCheckedParam<unsigned int>(
      "maximum_iterations", 80, "maximum_iterations>0", "Maximum bisection iterations.");
  params.addRangeCheckedParam<unsigned int>(
      "redistance_iterations",
      50,
      "redistance_iterations>0",
      "Number of second-order explicit TVD Runge-Kutta pseudo-time steps used to restore signed "
      "distance.");
  params.addParam<bool>("report", false, "Report the mass correction and interface shift.");
  params.addRelationshipManager("ElementSideNeighborLayers",
                                Moose::RelationshipManagerType::GEOMETRIC |
                                    Moose::RelationshipManagerType::ALGEBRAIC |
                                    Moose::RelationshipManagerType::COUPLING,
                                [](const InputParameters &, InputParameters & rm_params)
                                { rm_params.set<unsigned short>("layers") = 3; });
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_NONE);
  exec_enum = {EXEC_NONE};
  return params;
}

MassConservedLevelSetCorrector::MassConservedLevelSetCorrector(const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    BlockRestrictable(this),
    _system_name(
        UserObject::_subproblem.getVariable(0, getParam<VariableName>("level_set_variable"))
            .sys()
            .name()),
    _variable_name(getParam<VariableName>("level_set_variable")),
    _material_fraction_name(getParam<MooseFunctorName>("material_fraction")),
    _volumetric_face_flux_name(getParam<MooseFunctorName>("volumetric_face_flux")),
    _interface_width(getParam<Real>("interface_width")),
    _relative_tolerance(getParam<Real>("relative_tolerance")),
    _maximum_iterations(getParam<unsigned int>("maximum_iterations")),
    _redistance_iterations(getParam<unsigned int>("redistance_iterations")),
    _report(getParam<bool>("report")),
    _conserved_mass(
        declareRestartableData<Real>("conserved_mass", std::numeric_limits<Real>::quiet_NaN())),
    _level_set_variable(dynamic_cast<MooseLinearVariableFVReal &>(
        UserObject::_subproblem.getVariable(0, _variable_name))),
    _system(UserObject::_fe_problem.getLinearSystem(_level_set_variable.sys().number())),
    _system_number(_level_set_variable.sys().number()),
    _variable_number(_level_set_variable.number()),
    _material_density(getFunctor<Real>("material_density"))
{
}

void
MassConservedLevelSetCorrector::execute()
{
  correct();
}

void
MassConservedLevelSetCorrector::meshChanged()
{
  _redistance_stencils.clear();
  _redistance_value_cells.clear();
  _owned_redistance_cells.clear();
  _minimum_grid_spacing = 0.0;
}

const MassConservedLevelSetCorrector::CellStencil &
MassConservedLevelSetCorrector::stencil(const ElemInfo & elem_info) const
{
  const auto dof = elem_info.dofIndices()[_system_number][_variable_number];
  const auto stencil_it = _redistance_stencils.find(dof);
  if (stencil_it == _redistance_stencils.end())
    mooseError(name(), ": no redistance stencil is available for element ", elem_info.elem()->id());
  return stencil_it->second;
}

const MassConservedLevelSetCorrector::DirectionalStencil &
MassConservedLevelSetCorrector::direction(const ElemInfo & elem_info,
                                          const unsigned int axis,
                                          const bool positive) const
{
  return stencil(elem_info).directions[directionIndex(axis, positive)];
}

Real
MassConservedLevelSetCorrector::value(const ElemInfo & elem_info,
                                      const LevelSetValues & values) const
{
  const auto dof = elem_info.dofIndices()[_system_number][_variable_number];
  const auto value_it = values.find(dof);
  if (value_it == values.end())
    mooseError(name(), ": no redistance value is available for element ", elem_info.elem()->id());
  return value_it->second;
}

Real
MassConservedLevelSetCorrector::directionalValue(const ElemInfo & elem_info,
                                                 const unsigned int axis,
                                                 const bool positive,
                                                 const LevelSetValues & values) const
{
  const auto & cell_stencil = stencil(elem_info);
  const auto & directional_stencil = cell_stencil.directions[directionIndex(axis, positive)];
  if (directional_stencil.face_neighbors.empty())
    mooseError(name(),
               ": no directional redistance value is available at element ",
               elem_info.elem()->id(),
               ".");

  if (directional_stencil.conforming)
    return value(*directional_stencil.face_neighbors[0], values);

  const Real center_value = value(elem_info, values);
  Real reconstructed_value = center_value;
  Real minimum_value = center_value;
  Real maximum_value = center_value;
  mooseAssert(directional_stencil.reconstruction_weights.size() ==
                  cell_stencil.reconstruction_cells.size(),
              "The AMR reconstruction weights and patch cells must have equal sizes.");
  for (const auto sample : index_range(cell_stencil.reconstruction_cells))
  {
    const Real sample_value = value(*cell_stencil.reconstruction_cells[sample], values);
    reconstructed_value +=
        directional_stencil.reconstruction_weights[sample] * (sample_value - center_value);
    minimum_value = std::min(minimum_value, sample_value);
    maximum_value = std::max(maximum_value, sample_value);
  }

  // Limiting only the virtual coarse-fine value preserves the monotonicity required by the
  // Godunov Hamiltonian while retaining exact quadratic reconstruction away from extrema.
  return std::clamp(reconstructed_value, minimum_value, maximum_value);
}

void
MassConservedLevelSetCorrector::buildRedistanceStencils()
{
  if (!_redistance_stencils.empty())
    return;

  const unsigned int dimension = _fe_problem.mesh().dimension();
  if (dimension < 1 || dimension > 3)
    mooseError(name(),
               ": signed-distance redistancing requires a one-, two-, or three-dimensional "
               "Cartesian mesh.");

  const auto make_geometry_stencil = [this, dimension](const ElemInfo * const elem_info)
  {
    const auto * const elem = elem_info->elem();
    const dof_id_type dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (dof == DofObject::invalid_id)
      mooseError(name(),
                 ": the AMR redistance stencil has no level-set degree of "
                 "freedom on element ",
                 elem->id(),
                 ".");

    CellStencil cell_stencil;
    cell_stencil.elem_info = elem_info;
    std::array<bool, 2 * LIBMESH_DIM> direction_seen = {};
    auto inspect_face =
        [this, dimension, elem_info, &cell_stencil, &direction_seen](const Elem &,
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
                     ": AMR signed-distance redistancing requires axis-aligned Cartesian faces, "
                     "but the outward surface vector at element ",
                     elem_info->elem()->id(),
                     " is ",
                     surface_vector,
                     ".");

      const bool positive = surface_vector(axis) > 0.0;
      const auto direction_index = directionIndex(axis, positive);
      direction_seen[direction_index] = true;
      if (!neighbor || !hasBlocks(neighbor->subdomain_id()))
        return;

      if (std::abs(static_cast<int>(neighbor->level()) -
                   static_cast<int>(elem_info->elem()->level())) > 1)
        mooseError(name(),
                   ": redistancing requires a 2:1 balanced adaptive mesh, but elements ",
                   elem_info->elem()->id(),
                   " and ",
                   neighbor->id(),
                   " differ by more than one refinement level.");

      auto & directional_stencil = cell_stencil.directions[direction_index];
      directional_stencil.face_neighbors.push_back(&_fe_problem.mesh().elemInfo(neighbor->id()));
      directional_stencil.face_weights.push_back(surface_vector.norm());
    };

    const auto coord_type = _subproblem.getCoordSystem(elem->subdomain_id());
    Moose::FV::loopOverElemFaceInfo(*elem,
                                    _fe_problem.mesh(),
                                    inspect_face,
                                    coord_type,
                                    coord_type == Moose::COORD_RZ
                                        ? _subproblem.getAxisymmetricRadialCoord()
                                        : libMesh::invalid_uint);

    for (const auto axis : make_range(dimension))
      for (const bool positive : {false, true})
      {
        const auto direction_index = directionIndex(axis, positive);
        if (!direction_seen[direction_index])
          mooseError(name(),
                     ": element ",
                     elem->id(),
                     " has no Cartesian face in direction ",
                     positive ? "+" : "-",
                     axis,
                     ".");

        auto & directional_stencil = cell_stencil.directions[direction_index];
        if (directional_stencil.face_neighbors.empty())
          continue;

        Real weight_sum = 0.0;
        Point average_neighbor_centroid;
        for (const auto face : index_range(directional_stencil.face_neighbors))
        {
          const Real face_weight = directional_stencil.face_weights[face];
          weight_sum += face_weight;
          average_neighbor_centroid +=
              face_weight * directional_stencil.face_neighbors[face]->centroid();
        }
        average_neighbor_centroid /= weight_sum;
        directional_stencil.spacing =
            std::abs(average_neighbor_centroid(axis) - elem_info->centroid()(axis));
        if (directional_stencil.spacing <= 0.0)
          mooseError(name(), ": zero Cartesian neighbor spacing at element ", elem->id(), ".");

        directional_stencil.conforming =
            directional_stencil.face_neighbors.size() == 1 &&
            directional_stencil.face_neighbors[0]->elem()->level() == elem->level();
        if (directional_stencil.conforming)
        {
          const Point center_direction =
              directional_stencil.face_neighbors[0]->centroid() - elem_info->centroid();
          for (const auto component : make_range(dimension))
            if (component != axis &&
                std::abs(center_direction(component)) > 1e-10 * directional_stencil.spacing)
              directional_stencil.conforming = false;
        }
      }

    return std::make_pair(dof, cell_stencil);
  };

  std::vector<const ElemInfo *> stencil_cells;
  std::unordered_set<dof_id_type> stencil_dofs;
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    const auto * const elem = elem_info->elem();
    if (!hasBlocks(elem->subdomain_id()) || elem->processor_id() != processor_id())
      continue;

    _owned_redistance_cells.push_back(elem_info);
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    stencil_dofs.insert(dof);
    stencil_cells.push_back(elem_info);
  }

  // Neighbor curvatures are used by the ENO minmod correction, so geometry stencils are needed
  // for owned cells and their immediate face neighbors.
  for (std::size_t cell = 0; cell < stencil_cells.size(); ++cell)
  {
    const auto * const elem_info = stencil_cells[cell];
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (!_redistance_stencils.count(dof))
    {
      auto [stencil_dof, cell_stencil] = make_geometry_stencil(elem_info);
      mooseAssert(stencil_dof == dof, "Inconsistent redistance stencil degree of freedom.");
      for (const auto axis : make_range(dimension))
        for (const bool positive : {false, true})
          for (const auto * const neighbor :
               cell_stencil.directions[directionIndex(axis, positive)].face_neighbors)
          {
            const auto neighbor_dof = neighbor->dofIndices()[_system_number][_variable_number];
            if (stencil_dofs.insert(neighbor_dof).second)
              stencil_cells.push_back(neighbor);
          }
      _redistance_stencils.emplace(dof, std::move(cell_stencil));
    }

    // Stop after adding the first neighbor layer. Newly appended ghost cells receive geometry
    // stencils, but their neighbors are reconstruction values rather than additional stencils.
    if (cell + 1 == _owned_redistance_cells.size())
    {
      const auto first_neighbor_layer_size = stencil_cells.size();
      for (std::size_t ghost_cell = cell + 1; ghost_cell < first_neighbor_layer_size; ++ghost_cell)
      {
        const auto * const ghost_info = stencil_cells[ghost_cell];
        const auto ghost_dof = ghost_info->dofIndices()[_system_number][_variable_number];
        if (!_redistance_stencils.count(ghost_dof))
        {
          auto [built_dof, ghost_stencil] = make_geometry_stencil(ghost_info);
          mooseAssert(built_dof == ghost_dof,
                      "Inconsistent ghost redistance stencil degree of freedom.");
          _redistance_stencils.emplace(ghost_dof, std::move(ghost_stencil));
        }
      }
      break;
    }
  }

  const unsigned int coefficient_count = numberOfQuadraticCoefficients(dimension);
  for (auto & [dof, cell_stencil] : _redistance_stencils)
  {
    const auto * const center_info = cell_stencil.elem_info;
    std::vector<const ElemInfo *> frontier = {center_info};
    std::unordered_set<dof_id_type> patch_dofs = {dof};
    for (const auto layer : make_range(2))
    {
      (void)layer;
      std::vector<const ElemInfo *> next_frontier;
      for (const auto * const patch_cell : frontier)
      {
        const auto patch_dof = patch_cell->dofIndices()[_system_number][_variable_number];
        const auto stencil_it = _redistance_stencils.find(patch_dof);
        CellStencil temporary_stencil;
        const CellStencil * patch_stencil = nullptr;
        if (stencil_it != _redistance_stencils.end())
          patch_stencil = &stencil_it->second;
        else
        {
          temporary_stencil = make_geometry_stencil(patch_cell).second;
          patch_stencil = &temporary_stencil;
        }

        for (const auto axis : make_range(dimension))
          for (const bool positive : {false, true})
            for (const auto * const neighbor :
                 patch_stencil->directions[directionIndex(axis, positive)].face_neighbors)
            {
              const auto neighbor_dof = neighbor->dofIndices()[_system_number][_variable_number];
              if (patch_dofs.insert(neighbor_dof).second)
              {
                cell_stencil.reconstruction_cells.push_back(neighbor);
                next_frontier.push_back(neighbor);
              }
            }
      }
      frontier = std::move(next_frontier);
    }

    if (cell_stencil.reconstruction_cells.size() < coefficient_count)
      mooseError(name(),
                 ": the two-layer AMR patch around element ",
                 center_info->elem()->id(),
                 " contains only ",
                 cell_stencil.reconstruction_cells.size(),
                 " cells; ",
                 coefficient_count,
                 " are required for quadratic reconstruction.");

    const Real scale = center_info->elem()->hmin();
    const auto sample_count = cell_stencil.reconstruction_cells.size();
    DenseMatrix<Real> matrix(sample_count, coefficient_count);
    std::vector<Real> row_weights(sample_count);
    std::vector<std::vector<Real>> sample_bases(sample_count);
    for (const auto sample : index_range(cell_stencil.reconstruction_cells))
    {
      const Point offset =
          cell_stencil.reconstruction_cells[sample]->centroid() - center_info->centroid();
      sample_bases[sample] = quadraticBasis(offset, scale, dimension);
      row_weights[sample] = 1.0 / std::max(offset.norm() / scale, 1e-12);
      for (const auto coefficient : make_range(coefficient_count))
        matrix(sample, coefficient) = row_weights[sample] * sample_bases[sample][coefficient];
    }

    std::vector<std::vector<Real>> coefficient_weights(coefficient_count,
                                                       std::vector<Real>(sample_count));
    for (const auto sample : make_range(sample_count))
    {
      DenseVector<Real> rhs(sample_count);
      DenseVector<Real> solution(coefficient_count);
      rhs(sample) = row_weights[sample];
      matrix.svd_solve(rhs, solution);
      for (const auto coefficient : make_range(coefficient_count))
        coefficient_weights[coefficient][sample] = solution(coefficient);
    }

    for (const auto axis : make_range(dimension))
    {
      auto & weights = cell_stencil.second_derivative_weights[axis];
      weights.resize(sample_count);
      for (const auto sample : make_range(sample_count))
        weights[sample] = coefficient_weights[dimension + axis][sample] / (scale * scale);
    }

    for (const auto axis : make_range(dimension))
      for (const bool positive : {false, true})
      {
        auto & directional_stencil = cell_stencil.directions[directionIndex(axis, positive)];
        if (directional_stencil.face_neighbors.empty() || directional_stencil.conforming)
          continue;

        Point target_offset;
        target_offset(axis) = (positive ? 1.0 : -1.0) * directional_stencil.spacing;
        const auto target_basis = quadraticBasis(target_offset, scale, dimension);
        directional_stencil.reconstruction_weights.resize(sample_count);
        for (const auto sample : make_range(sample_count))
          for (const auto coefficient : make_range(coefficient_count))
            directional_stencil.reconstruction_weights[sample] +=
                target_basis[coefficient] * coefficient_weights[coefficient][sample];

        for (const auto coefficient : make_range(coefficient_count))
        {
          Real reproduced_basis = 0.0;
          for (const auto sample : make_range(sample_count))
            reproduced_basis += directional_stencil.reconstruction_weights[sample] *
                                sample_bases[sample][coefficient];
          if (std::abs(reproduced_basis - target_basis[coefficient]) > 1e-8)
            mooseError(name(),
                       ": singular quadratic AMR reconstruction at element ",
                       center_info->elem()->id(),
                       " in coefficient ",
                       coefficient,
                       ".");
        }
      }

    _redistance_value_cells.emplace(dof, center_info);
    for (const auto * const sample : cell_stencil.reconstruction_cells)
    {
      const auto sample_dof = sample->dofIndices()[_system_number][_variable_number];
      _redistance_value_cells.emplace(sample_dof, sample);
    }
  }

  Real local_minimum_spacing = std::numeric_limits<Real>::max();
  for (const auto * const elem_info : _owned_redistance_cells)
    for (const auto axis : make_range(dimension))
      for (const bool positive : {false, true})
      {
        const Real spacing = direction(*elem_info, axis, positive).spacing;
        if (spacing > 0.0)
          local_minimum_spacing = std::min(local_minimum_spacing, spacing);
      }
  _communicator.min(local_minimum_spacing);
  if (!std::isfinite(local_minimum_spacing))
    mooseError(name(), ": no internal Cartesian spacing is available for redistancing.");
  _minimum_grid_spacing = local_minimum_spacing;
}

Real
MassConservedLevelSetCorrector::secondDerivative(const ElemInfo & elem_info,
                                                 const unsigned int axis,
                                                 const LevelSetValues & values) const
{
  const auto & cell_stencil = stencil(elem_info);
  const auto & negative = direction(elem_info, axis, false);
  const auto & positive = direction(elem_info, axis, true);
  const Real center_value = value(elem_info, values);
  if (!negative.face_neighbors.empty() && !positive.face_neighbors.empty())
  {
    const Real negative_slope =
        (center_value - directionalValue(elem_info, axis, false, values)) / negative.spacing;
    const Real positive_slope =
        (directionalValue(elem_info, axis, true, values) - center_value) / positive.spacing;
    return 2.0 * (positive_slope - negative_slope) / (negative.spacing + positive.spacing);
  }

  Real reconstructed_derivative = 0.0;
  const auto & weights = cell_stencil.second_derivative_weights[axis];
  for (const auto sample : index_range(cell_stencil.reconstruction_cells))
    reconstructed_derivative +=
        weights[sample] *
        (value(*cell_stencil.reconstruction_cells[sample], values) - center_value);
  return reconstructed_derivative;
}

Real
MassConservedLevelSetCorrector::subcellDistance(const ElemInfo & elem_info,
                                                const unsigned int axis,
                                                const bool positive,
                                                const LevelSetValues & initial_values) const
{
  const auto & directional_stencil = direction(elem_info, axis, positive);
  mooseAssert(!directional_stencil.face_neighbors.empty(),
              "A subcell interface distance requires an adjacent cell.");
  const Real center_value = value(elem_info, initial_values);
  const Real adjacent_value = directionalValue(elem_info, axis, positive, initial_values);
  const Real spacing = directional_stencil.spacing;
  const Real linear_distance = spacing * std::abs(center_value / (center_value - adjacent_value));

  Real neighbor_curvature = 0.0;
  Real neighbor_weight = 0.0;
  for (const auto face : index_range(directional_stencil.face_neighbors))
  {
    neighbor_curvature +=
        directional_stencil.face_weights[face] *
        secondDerivative(*directional_stencil.face_neighbors[face], axis, initial_values);
    neighbor_weight += directional_stencil.face_weights[face];
  }
  neighbor_curvature /= neighbor_weight;
  const Real eno_curvature =
      minmod(secondDerivative(elem_info, axis, initial_values), neighbor_curvature);
  const Real curvature_tolerance =
      1e-12 * std::max({1.0, std::abs(center_value), std::abs(adjacent_value)}) /
      (spacing * spacing);
  if (std::abs(eno_curvature) <= curvature_tolerance)
    return linear_distance;

  const Real slope =
      (adjacent_value - center_value - 0.5 * eno_curvature * spacing * spacing) / spacing;
  const Real discriminant = slope * slope - 2.0 * eno_curvature * center_value;
  if (discriminant < 0.0)
    return linear_distance;

  const Real square_root = std::sqrt(discriminant);
  const std::array<Real, 2> roots = {(-slope - square_root) / eno_curvature,
                                     (-slope + square_root) / eno_curvature};
  Real quadratic_distance = std::numeric_limits<Real>::max();
  for (const Real root : roots)
    if (root > 0.0 && root <= spacing)
      quadratic_distance = std::min(quadratic_distance, root);
  return std::isfinite(quadratic_distance) ? quadratic_distance : linear_distance;
}

Real
MassConservedLevelSetCorrector::oneSidedDerivative(const ElemInfo & elem_info,
                                                   const unsigned int axis,
                                                   const bool positive,
                                                   const LevelSetValues & values,
                                                   const LevelSetValues & initial_values) const
{
  const auto & directional_stencil = direction(elem_info, axis, positive);
  mooseAssert(!directional_stencil.face_neighbors.empty(),
              "A one-sided derivative requires an adjacent cell.");
  const Real center_value = value(elem_info, values);
  const Real adjacent_value = directionalValue(elem_info, axis, positive, values);
  const Real center_initial_value = value(elem_info, initial_values);
  const Real adjacent_initial_value = directionalValue(elem_info, axis, positive, initial_values);
  const bool crosses_interface =
      center_initial_value != 0.0 && center_initial_value * adjacent_initial_value <= 0.0;
  const Real distance = crosses_interface
                            ? subcellDistance(elem_info, axis, positive, initial_values)
                            : directional_stencil.spacing;
  Real neighbor_curvature = 0.0;
  Real neighbor_weight = 0.0;
  for (const auto face : index_range(directional_stencil.face_neighbors))
  {
    neighbor_curvature += directional_stencil.face_weights[face] *
                          secondDerivative(*directional_stencil.face_neighbors[face], axis, values);
    neighbor_weight += directional_stencil.face_weights[face];
  }
  neighbor_curvature /= neighbor_weight;
  const Real eno_curvature = minmod(secondDerivative(elem_info, axis, values), neighbor_curvature);

  if (positive)
    return (crosses_interface ? -center_value : adjacent_value - center_value) / distance -
           0.5 * distance * eno_curvature;

  return (crosses_interface ? center_value : center_value - adjacent_value) / distance +
         0.5 * distance * eno_curvature;
}

Real
MassConservedLevelSetCorrector::godunovHamiltonian(const ElemInfo & elem_info,
                                                   const Real initial_value,
                                                   const LevelSetValues & values,
                                                   const LevelSetValues & initial_values) const
{
  Real squared_hamiltonian = 0.0;
  for (const auto axis : make_range(_fe_problem.mesh().dimension()))
  {
    const bool has_negative = !direction(elem_info, axis, false).face_neighbors.empty();
    const bool has_positive = !direction(elem_info, axis, true).face_neighbors.empty();
    if (!has_negative && !has_positive)
      mooseError(name(),
                 ": no one-sided redistance derivative is available at element ",
                 elem_info.elem()->id(),
                 " in direction ",
                 axis,
                 ".");

    Real positive_derivative =
        has_positive ? oneSidedDerivative(elem_info, axis, true, values, initial_values) : 0.0;
    Real negative_derivative =
        has_negative ? oneSidedDerivative(elem_info, axis, false, values, initial_values) : 0.0;
    if (!has_positive)
      positive_derivative = negative_derivative;
    if (!has_negative)
      negative_derivative = positive_derivative;

    const Real selected = initial_value >= 0.0
                              ? std::max(Utility::pow<2>(std::min(positive_derivative, 0.0)),
                                         Utility::pow<2>(std::max(negative_derivative, 0.0)))
                              : std::max(Utility::pow<2>(std::max(positive_derivative, 0.0)),
                                         Utility::pow<2>(std::min(negative_derivative, 0.0)));
    squared_hamiltonian += selected;
  }
  return std::sqrt(squared_hamiltonian);
}

Real
MassConservedLevelSetCorrector::localPseudoTimeStep(const ElemInfo & elem_info,
                                                    const LevelSetValues & initial_values) const
{
  const unsigned int dimension = _fe_problem.mesh().dimension();
  Real minimum_distance = std::numeric_limits<Real>::max();
  const Real center_value = value(elem_info, initial_values);
  for (const auto axis : make_range(dimension))
    for (const bool positive : {false, true})
    {
      const auto & directional_stencil = direction(elem_info, axis, positive);
      if (!directional_stencil.face_neighbors.empty())
      {
        const Real adjacent_value = directionalValue(elem_info, axis, positive, initial_values);
        const Real distance = center_value != 0.0 && center_value * adjacent_value <= 0.0
                                  ? subcellDistance(elem_info, axis, positive, initial_values)
                                  : directional_stencil.spacing;
        minimum_distance = std::min(minimum_distance, distance);
      }
    }

  return minimum_distance / dimension;
}

void
MassConservedLevelSetCorrector::synchronizeLevelSet(LevelSetValues & values,
                                                    const Moose::StateArg & state)
{
  auto & algebraic_solution = _system.solutionState(state.state, state.iteration_type);
  auto & current_local_solution = *_system.system().current_local_solution;

  for (const auto * const elem_info : _owned_redistance_cells)
  {
    const dof_id_type dof = elem_info->dofIndices()[_system_number][_variable_number];
    algebraic_solution.set(dof, value(*elem_info, values));
  }
  algebraic_solution.close();
  algebraic_solution.localize(current_local_solution, _system.dofMap().get_send_list());
  if (state.state == 0)
    _system.setSolution(current_local_solution);

  for (const auto & [dof, elem_info] : _redistance_value_cells)
    values[dof] = current_local_solution(dof);
}

Real
MassConservedLevelSetCorrector::redistance(const Moose::StateArg & state)
{
  buildRedistanceStencils();

  LevelSetValues initial_values;
  LevelSetValues values;
  for (const auto & [dof, elem_info] : _redistance_value_cells)
  {
    const Real initial_value = _level_set_variable.getElemValue(*elem_info, state);
    initial_values.emplace(dof, initial_value);
    values.emplace(dof, initial_value);
  }

  const unsigned int dimension = _fe_problem.mesh().dimension();
  std::unordered_set<dof_id_type> active_cells;
  std::unordered_set<dof_id_type> convergence_cells;
  for (const auto * const elem_info : _owned_redistance_cells)
  {
    const Real center_value = value(*elem_info, initial_values);
    Real squared_gradient = 0.0;
    Real local_maximum_spacing = 0.0;
    bool touches_interface = center_value == 0.0;
    for (const auto axis : make_range(dimension))
    {
      const auto & negative = direction(*elem_info, axis, false);
      const auto & positive = direction(*elem_info, axis, true);
      for (const auto * const neighbor : negative.face_neighbors)
        touches_interface =
            touches_interface || center_value * value(*neighbor, initial_values) <= 0.0;
      for (const auto * const neighbor : positive.face_neighbors)
        touches_interface =
            touches_interface || center_value * value(*neighbor, initial_values) <= 0.0;

      Real derivative = 0.0;
      if (!negative.face_neighbors.empty() && !positive.face_neighbors.empty())
      {
        const Real negative_value = directionalValue(*elem_info, axis, false, initial_values);
        const Real positive_value = directionalValue(*elem_info, axis, true, initial_values);
        derivative = -positive.spacing /
                         (negative.spacing * (negative.spacing + positive.spacing)) *
                         negative_value +
                     (positive.spacing - negative.spacing) / (negative.spacing * positive.spacing) *
                         center_value +
                     negative.spacing / (positive.spacing * (negative.spacing + positive.spacing)) *
                         positive_value;
      }
      else if (!positive.face_neighbors.empty())
        derivative = (directionalValue(*elem_info, axis, true, initial_values) - center_value) /
                     positive.spacing;
      else if (!negative.face_neighbors.empty())
        derivative = (center_value - directionalValue(*elem_info, axis, false, initial_values)) /
                     negative.spacing;
      squared_gradient += derivative * derivative;
      local_maximum_spacing = std::max({local_maximum_spacing, negative.spacing, positive.spacing});
    }

    const Real gradient_magnitude = std::sqrt(squared_gradient);
    const Real estimated_distance = gradient_magnitude > std::numeric_limits<Real>::epsilon()
                                        ? std::abs(center_value) / gradient_magnitude
                                        : std::numeric_limits<Real>::max();
    const Real convergence_band_width = _interface_width + 2.0 * local_maximum_spacing;
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    // CASL advances every node during each redistance stage. Restricting the update to a narrow
    // band introduces a frozen internal boundary whose values enter the second-order stencil and
    // contaminate the interface under repeated RK2 sweeps.
    active_cells.insert(dof);
    if (touches_interface || estimated_distance <= convergence_band_width)
      convergence_cells.insert(dof);
  }

  std::unordered_map<dof_id_type, Real> pseudo_time_steps;
  for (const auto * const elem_info : _owned_redistance_cells)
  {
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (active_cells.count(dof))
      pseudo_time_steps.emplace(dof, localPseudoTimeStep(*elem_info, initial_values));
  }

  const auto reinitialization_operator =
      [this, &initial_values](const ElemInfo & elem_info, const LevelSetValues & stage_values)
  {
    const Real initial_value = value(elem_info, initial_values);
    return strictSign(initial_value) *
           (godunovHamiltonian(elem_info, initial_value, stage_values, initial_values) - 1.0);
  };

  const auto forward_euler_step =
      [this,
       &active_cells,
       &initial_values,
       &pseudo_time_steps,
       &reinitialization_operator,
       &state](const LevelSetValues & current_values, LevelSetValues & next_values)
  {
    next_values = current_values;
    for (const auto * const elem_info : _owned_redistance_cells)
    {
      const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
      if (!active_cells.count(dof))
        continue;

      const Real initial_value = value(*elem_info, initial_values);
      if (initial_value == 0.0)
      {
        next_values[dof] = 0.0;
        continue;
      }

      Real new_value =
          value(*elem_info, current_values) -
          pseudo_time_steps.at(dof) * reinitialization_operator(*elem_info, current_values);
      if (initial_value * new_value < 0.0)
        new_value *= -1.0;
      next_values[dof] = new_value;
    }

    synchronizeLevelSet(next_values, state);
  };

  Real final_normalized_update = 0.0;
  for (const auto pseudo_step : make_range(_redistance_iterations))
  {
    (void)pseudo_step;
    const LevelSetValues previous_values(values);
    LevelSetValues first_stage;
    forward_euler_step(previous_values, first_stage);
    LevelSetValues second_stage;
    forward_euler_step(first_stage, second_stage);

    for (const auto * const elem_info : _owned_redistance_cells)
    {
      const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
      if (active_cells.count(dof))
        values[dof] = 0.5 * (value(*elem_info, previous_values) + value(*elem_info, second_stage));
    }
    synchronizeLevelSet(values, state);

    Real local_maximum_update = 0.0;
    for (const auto * const elem_info : _owned_redistance_cells)
    {
      const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
      if (!convergence_cells.count(dof))
        continue;
      local_maximum_update =
          std::max(local_maximum_update,
                   std::abs(value(*elem_info, values) - value(*elem_info, previous_values)));
    }
    _communicator.max(local_maximum_update);
    final_normalized_update = local_maximum_update / _minimum_grid_spacing;
  }

  unsigned int changed_phase_signs = 0;
  for (const auto * const elem_info : _owned_redistance_cells)
  {
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (active_cells.count(dof) &&
        value(*elem_info, values) * value(*elem_info, initial_values) < 0.0)
      ++changed_phase_signs;
  }
  _communicator.sum(changed_phase_signs);
  if (changed_phase_signs)
    mooseError(name(),
               ": ENO redistancing changed the phase sign in ",
               changed_phase_signs,
               " narrow-band cells.");

  if (_report)
    _console << name() << ": TVD RK2 redistance completed " << _redistance_iterations
             << " pseudo-time iterations over " << convergence_cells.size()
             << " owned core cells and " << active_cells.size()
             << " owned active cells with normalized update " << final_normalized_update
             << std::endl;
  return final_normalized_update;
}

Real
MassConservedLevelSetCorrector::regularizedHeaviside(const Real phi) const
{
  if (phi <= -_interface_width)
    return 0.0;
  if (phi >= _interface_width)
    return 1.0;

  const Real normalized_phi = phi / _interface_width;
  return 0.5 * (1.0 + normalized_phi + std::sin(libMesh::pi * normalized_phi) / libMesh::pi);
}

Real
MassConservedLevelSetCorrector::materialMass(const Real shift, const Moose::StateArg & state) const
{
  Real mass = 0.0;
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    const auto * const elem = elem_info->elem();
    if (!hasBlocks(elem->subdomain_id()) || elem->processor_id() != processor_id())
      continue;

    const Real phi = _level_set_variable.getElemValue(*elem_info, state);
    const Real density = _material_density(makeElemArg(elem), state);
    if (!std::isfinite(phi) || !std::isfinite(density) || density <= 0.0)
      mooseError(name(),
                 ": nonphysical level-set material state in element ",
                 elem->id(),
                 ": phi=",
                 phi,
                 ", density=",
                 density,
                 ".");

    mass += regularizedHeaviside(phi + shift) * density * elem_info->volume() *
            elem_info->coordFactor();
  }

  _communicator.sum(mass);
  return mass;
}

Real
MassConservedLevelSetCorrector::massConservingShift(const Real target_mass,
                                                    const Real uncorrected_mass) const
{
  const Real mass_tolerance =
      _relative_tolerance * std::max(target_mass, std::numeric_limits<Real>::min());

  if (std::abs(uncorrected_mass - target_mass) <= mass_tolerance)
    return 0.0;

  Real local_extent = _interface_width;
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
    if (hasBlocks(elem_info->subdomain_id()))
      local_extent =
          std::max(local_extent,
                   std::abs(_level_set_variable.getElemValue(*elem_info, Moose::currentState())) +
                       2.0 * _interface_width);
  _communicator.max(local_extent);

  Real lower_shift = -local_extent;
  Real upper_shift = local_extent;
  const Real lower_mass = materialMass(lower_shift, Moose::currentState());
  const Real upper_mass = materialMass(upper_shift, Moose::currentState());
  if (target_mass < lower_mass - mass_tolerance || target_mass > upper_mass + mass_tolerance)
    mooseError(name(),
               ": conserved material mass ",
               target_mass,
               " is outside the current level-set representation range [",
               lower_mass,
               ", ",
               upper_mass,
               "].");

  Real shift = 0.0;
  Real corrected_mass = uncorrected_mass;
  for (const auto iteration : make_range(_maximum_iterations))
  {
    (void)iteration;
    shift = 0.5 * (lower_shift + upper_shift);
    corrected_mass = materialMass(shift, Moose::currentState());
    if (std::abs(corrected_mass - target_mass) <= mass_tolerance)
      break;

    if (corrected_mass < target_mass)
      lower_shift = shift;
    else
      upper_shift = shift;
  }

  if (std::abs(corrected_mass - target_mass) > mass_tolerance)
    mooseError(name(),
               ": failed to satisfy the material-mass constraint. Target mass: ",
               target_mass,
               "; corrected mass: ",
               corrected_mass,
               "; interface shift: ",
               shift,
               ".");

  return shift;
}

Real
MassConservedLevelSetCorrector::correctMassOnly()
{
  if (!std::isfinite(_conserved_mass))
    _conserved_mass = materialMass(0.0, Moose::oldState());

  const Real target_mass = _conserved_mass;
  const Real uncorrected_mass = materialMass(0.0, Moose::currentState());
  const Real shift = massConservingShift(target_mass, uncorrected_mass);
  if (shift != 0.0)
    applyShift(shift);
  else if (!_correcting_after_mesh_change)
    _system.computeGradients();

  if (_report)
    _console << name() << ": level-set mass " << uncorrected_mass << " -> "
             << materialMass(0.0, Moose::currentState()) << " (target " << target_mass
             << "), uniform shift " << shift << std::endl;
  return shift;
}

Real
MassConservedLevelSetCorrector::correct()
{
  LevelSetValues pre_correction_values;
  buildRedistanceStencils();
  for (const auto * const elem_info : _owned_redistance_cells)
  {
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    pre_correction_values.emplace(
        dof, _level_set_variable.getElemValue(*elem_info, Moose::currentState()));
  }

  redistance();
  const Real shift = correctMassOnly();
  _last_normalized_update = 0.0;
  for (const auto * const elem_info : _owned_redistance_cells)
  {
    const auto dof = elem_info->dofIndices()[_system_number][_variable_number];
    _last_normalized_update =
        std::max(_last_normalized_update,
                 std::abs(_level_set_variable.getElemValue(*elem_info, Moose::currentState()) -
                          pre_correction_values.at(dof)) /
                     _minimum_grid_spacing);
  }
  _communicator.max(_last_normalized_update);
  if (_report)
    _console << name() << ": normalized total redistance/mass correction "
             << _last_normalized_update << std::endl;
  return shift;
}

void
MassConservedLevelSetCorrector::correctAfterMeshChange()
{
  meshChanged();
  _correcting_after_mesh_change = true;

  unsigned int number_of_states = 0;
  while (_system.hasSolutionState(number_of_states))
    ++number_of_states;

  // Correct old states before current because arbitrary-state synchronization uses the current
  // local vector as ghosted scratch storage. The current state is restored before its correction.
  for (unsigned int state_index = number_of_states; state_index-- > 1;)
  {
    const Moose::StateArg state{state_index, Moose::SolutionIterationType::Time};
    const Real target_mass = materialMass(0.0, state);
    redistance(state);
    const Real uncorrected_mass = materialMass(0.0, state);
    const Real shift = massConservingShift(target_mass, uncorrected_mass);
    if (shift != 0.0)
      applyShift(shift, state);

    if (_report)
      _console << name() << ": post-remesh state " << state_index << " level-set mass "
               << uncorrected_mass << " -> " << materialMass(0.0, state) << " (target "
               << target_mass << "), uniform shift " << shift << std::endl;
  }

  restoreCurrentSolution();
  correct();
  _correcting_after_mesh_change = false;
}

void
MassConservedLevelSetCorrector::restoreCurrentSolution()
{
  auto & current_local_solution = *_system.system().current_local_solution;
  _system.solution().localize(current_local_solution, _system.dofMap().get_send_list());
  _system.setSolution(current_local_solution);
  if (!_correcting_after_mesh_change)
    _system.computeGradients();
}

void
MassConservedLevelSetCorrector::applyShift(const Real shift, const Moose::StateArg & state)
{
  auto & algebraic_solution = _system.solutionState(state.state, state.iteration_type);

  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    const auto * const elem = elem_info->elem();
    if (!hasBlocks(elem->subdomain_id()) || elem->processor_id() != processor_id())
      continue;

    const dof_id_type dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (dof == DofObject::invalid_id || dof < algebraic_solution.first_local_index() ||
        dof >= algebraic_solution.last_local_index())
      continue;

    const Real shifted_phi = _level_set_variable.getElemValue(*elem_info, state) + shift;
    algebraic_solution.set(dof, shifted_phi);
  }

  algebraic_solution.close();
  if (state.state == 0)
    restoreCurrentSolution();
}
