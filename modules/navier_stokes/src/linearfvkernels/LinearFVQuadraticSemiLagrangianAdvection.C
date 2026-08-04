//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVQuadraticSemiLagrangianAdvection.h"

#include "FEProblemBase.h"
#include "FVUtils.h"
#include "MooseMesh.h"
#include "SubProblem.h"

#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/remote_elem.h"

#include <algorithm>
#include <cmath>
#include <unordered_set>

registerMooseObject("NavierStokesApp", LinearFVQuadraticSemiLagrangianAdvection);

InputParameters
LinearFVQuadraticSemiLagrangianAdvection::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Advects a linear finite-volume scalar by tracing a midpoint characteristic backward and "
      "quadratically reconstructing the old solution at its departure point.");
  params.addRequiredParam<MooseFunctorName>("u", "Velocity component in the x direction.");
  params.addParam<MooseFunctorName>("v", "0", "Velocity component in the y direction.");
  params.addParam<MooseFunctorName>("w", "0", "Velocity component in the z direction.");
  params.addRangeCheckedParam<Real>(
      "maximum_courant_number",
      1.0,
      "maximum_courant_number>0 & maximum_courant_number<=8",
      "Maximum permitted characteristic displacement divided by the arrival-cell minimum size.");
  params.set<unsigned short>("ghost_layers") = 3;
  params.set<bool>("use_point_neighbors") = true;
  return params;
}

LinearFVQuadraticSemiLagrangianAdvection::LinearFVQuadraticSemiLagrangianAdvection(
    const InputParameters & params)
  : LinearFVElementalKernel(params),
    _dimension(_fe_problem.mesh().dimension()),
    _velocity{{&getFunctor<Real>("u"), &getFunctor<Real>("v"), &getFunctor<Real>("w")}},
    _maximum_courant_number(getParam<Real>("maximum_courant_number"))
{
  if (_dimension < 1 || _dimension > 3)
    mooseError(name(), ": semi-Lagrangian advection requires a 1D, 2D, or 3D mesh.");

  const auto required_ghost_layers =
      static_cast<unsigned short>(std::ceil(_maximum_courant_number)) + 2;
  if (getParam<unsigned short>("ghost_layers") < required_ghost_layers)
    paramError("ghost_layers",
               "At least ",
               required_ghost_layers,
               " point-neighbor layers are required for the requested maximum Courant number "
               "and the two-layer quadratic reconstruction.");
  if (!getParam<bool>("use_point_neighbors"))
    paramError("use_point_neighbors",
               "Point-neighbor ghosting is required because a characteristic may cross a cell "
               "corner.");
}

void
LinearFVQuadraticSemiLagrangianAdvection::initialSetup()
{
  resetGeometry();
}

void
LinearFVQuadraticSemiLagrangianAdvection::meshChanged()
{
  resetGeometry();
}

void
LinearFVQuadraticSemiLagrangianAdvection::resetGeometry()
{
  _reconstructions.clear();
  _point_locator = _fe_problem.mesh().getPointLocator();
  _point_locator->enable_out_of_mesh_mode();
  for (const auto axis : make_range(_dimension))
  {
    _domain_min[axis] = _fe_problem.mesh().getMinInDimension(axis);
    _domain_max[axis] = _fe_problem.mesh().getMaxInDimension(axis);
  }
}

Real
LinearFVQuadraticSemiLagrangianAdvection::computeMatrixContribution()
{
  return _current_elem_volume;
}

Real
LinearFVQuadraticSemiLagrangianAdvection::computeRightHandSideContribution()
{
  return _current_elem_volume *
         reconstructedVariable(departurePoint(*_current_elem_info), Moose::oldState());
}

const ElemInfo &
LinearFVQuadraticSemiLagrangianAdvection::containingElement(const Point & point) const
{
  const Elem * elem =
      blockRestricted() ? (*_point_locator)(point, &blockIDs()) : (*_point_locator)(point);
  if (!elem || elem == libMesh::remote_elem)
    mooseError(name(),
               ": departure point ",
               point,
               " is not in the local or ghosted mesh. Increase ghost_layers or reduce the "
               "maximum Courant number.");
  return _fe_problem.mesh().elemInfo(elem->id());
}

std::vector<const ElemInfo *>
LinearFVQuadraticSemiLagrangianAdvection::faceNeighbors(const ElemInfo & center) const
{
  std::vector<const ElemInfo *> neighbors;
  const Elem * const elem = center.elem();
  auto inspect_face = [this, &neighbors](const Elem &,
                                         const Elem * const neighbor,
                                         const FaceInfo *,
                                         const Point &,
                                         const Real,
                                         const bool)
  {
    if (neighbor && neighbor != libMesh::remote_elem && hasBlocks(neighbor->subdomain_id()))
      neighbors.push_back(&_fe_problem.mesh().elemInfo(neighbor->id()));
  };

  const auto coord_type = _subproblem.getCoordSystem(elem->subdomain_id());
  Moose::FV::loopOverElemFaceInfo(*elem,
                                  _fe_problem.mesh(),
                                  inspect_face,
                                  coord_type,
                                  coord_type == Moose::COORD_RZ
                                      ? _subproblem.getAxisymmetricRadialCoord()
                                      : libMesh::invalid_uint);
  return neighbors;
}

std::vector<Real>
LinearFVQuadraticSemiLagrangianAdvection::quadraticBasis(const Point & offset,
                                                         const Real scale) const
{
  const unsigned int coefficient_count = _dimension * (_dimension + 3) / 2;
  std::vector<Real> basis(coefficient_count);
  std::array<Real, LIBMESH_DIM> coordinate = {};
  for (const auto axis : make_range(_dimension))
  {
    coordinate[axis] = offset(axis) / scale;
    basis[axis] = coordinate[axis];
    basis[_dimension + axis] = 0.5 * coordinate[axis] * coordinate[axis];
  }

  unsigned int coefficient = 2 * _dimension;
  for (const auto first : make_range(_dimension))
    for (const auto second : make_range(first + 1, _dimension))
      basis[coefficient++] = coordinate[first] * coordinate[second];
  return basis;
}

const LinearFVQuadraticSemiLagrangianAdvection::Reconstruction &
LinearFVQuadraticSemiLagrangianAdvection::reconstruction(const ElemInfo & center)
{
  const auto center_id = center.elem()->id();
  if (const auto it = _reconstructions.find(center_id); it != _reconstructions.end())
    return it->second;

  Reconstruction result;
  result.scale = center.elem()->hmin();

  std::vector<const ElemInfo *> frontier = {&center};
  std::unordered_set<dof_id_type> patch_ids = {center_id};
  for (const auto layer : make_range(2))
  {
    (void)layer;
    std::vector<const ElemInfo *> next_frontier;
    for (const auto * const patch_cell : frontier)
      for (const auto * const neighbor : faceNeighbors(*patch_cell))
        if (patch_ids.insert(neighbor->elem()->id()).second)
        {
          result.samples.push_back(neighbor);
          next_frontier.push_back(neighbor);
        }
    frontier = std::move(next_frontier);
  }

  const unsigned int coefficient_count = _dimension * (_dimension + 3) / 2;
  if (result.samples.size() < coefficient_count)
    mooseError(name(),
               ": the two-layer patch around element ",
               center_id,
               " contains ",
               result.samples.size(),
               " cells, but ",
               coefficient_count,
               " are required for quadratic reconstruction.");

  const auto sample_count = result.samples.size();
  DenseMatrix<Real> matrix(sample_count, coefficient_count);
  std::vector<Real> row_weights(sample_count);
  for (const auto sample : index_range(result.samples))
  {
    const Point offset = result.samples[sample]->centroid() - center.centroid();
    const auto basis = quadraticBasis(offset, result.scale);
    row_weights[sample] = 1.0 / std::max(offset.norm() / result.scale, 1e-12);
    for (const auto coefficient : make_range(coefficient_count))
      matrix(sample, coefficient) = row_weights[sample] * basis[coefficient];
  }

  result.coefficient_weights.assign(coefficient_count, std::vector<Real>(sample_count));
  for (const auto sample : make_range(sample_count))
  {
    DenseVector<Real> rhs(sample_count);
    DenseVector<Real> solution(coefficient_count);
    rhs(sample) = row_weights[sample];
    matrix.svd_solve(rhs, solution);
    for (const auto coefficient : make_range(coefficient_count))
      result.coefficient_weights[coefficient][sample] = solution(coefficient);
  }

  return _reconstructions.emplace(center_id, std::move(result)).first->second;
}

Real
LinearFVQuadraticSemiLagrangianAdvection::reconstructedVariable(const Point & point,
                                                                const Moose::StateArg & state)
{
  const auto & center = containingElement(point);
  const auto & patch = reconstruction(center);
  const Real center_value = _var.getElemValue(center, state);
  const auto target_basis = quadraticBasis(point - center.centroid(), patch.scale);
  Real reconstructed_value = center_value;
  for (const auto sample : index_range(patch.samples))
  {
    Real sample_weight = 0.0;
    for (const auto coefficient : index_range(target_basis))
      sample_weight += target_basis[coefficient] * patch.coefficient_weights[coefficient][sample];
    reconstructed_value +=
        sample_weight * (_var.getElemValue(*patch.samples[sample], state) - center_value);
  }
  return reconstructed_value;
}

Real
LinearFVQuadraticSemiLagrangianAdvection::reconstructedFunctor(const Moose::Functor<Real> & functor,
                                                               const Point & point,
                                                               const Moose::StateArg & state)
{
  const auto & center = containingElement(point);
  const auto & patch = reconstruction(center);
  const Real center_value = functor(makeElemArg(center.elem()), state);
  const auto target_basis = quadraticBasis(point - center.centroid(), patch.scale);
  Real reconstructed_value = center_value;
  for (const auto sample : index_range(patch.samples))
  {
    Real sample_weight = 0.0;
    for (const auto coefficient : index_range(target_basis))
      sample_weight += target_basis[coefficient] * patch.coefficient_weights[coefficient][sample];
    reconstructed_value +=
        sample_weight * (functor(makeElemArg(patch.samples[sample]->elem()), state) - center_value);
  }
  return reconstructed_value;
}

Point
LinearFVQuadraticSemiLagrangianAdvection::boundedPoint(Point point) const
{
  for (const auto axis : make_range(_dimension))
    if (_fe_problem.mesh().isTranslatedPeriodic(_var, axis))
    {
      const Real length = _domain_max[axis] - _domain_min[axis];
      point(axis) = _domain_min[axis] +
                    std::fmod(std::fmod(point(axis) - _domain_min[axis], length) + length, length);
    }
    else
      point(axis) = std::clamp(point(axis), _domain_min[axis], _domain_max[axis]);
  return point;
}

Point
LinearFVQuadraticSemiLagrangianAdvection::departurePoint(const ElemInfo & arrival)
{
  const auto state = Moose::oldState();
  const Point arrival_point = arrival.centroid();
  Point midpoint = arrival_point;
  for (const auto axis : make_range(_dimension))
    midpoint(axis) -= 0.5 * _dt * reconstructedFunctor(*_velocity[axis], arrival_point, state);
  midpoint = boundedPoint(midpoint);

  Point departure = arrival_point;
  for (const auto axis : make_range(_dimension))
    departure(axis) -= _dt * reconstructedFunctor(*_velocity[axis], midpoint, state);
  departure = boundedPoint(departure);

  const Real displacement = (departure - arrival_point).norm();
  const Real permitted_displacement = _maximum_courant_number * arrival.elem()->hmin();
  if (displacement > permitted_displacement * (1.0 + 1e-8))
    mooseError(name(),
               ": characteristic displacement ",
               displacement,
               " exceeds maximum_courant_number * hmin = ",
               permitted_displacement,
               " on element ",
               arrival.elem()->id(),
               ". Reduce dt or increase maximum_courant_number and ghost_layers.");
  return departure;
}
