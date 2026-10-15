//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosGreenGaussGradient.h"

#include "KokkosFVMesh.h"
#include "KokkosLinearFVBoundaryCondition.h"
#include "KokkosSystem.h"

#include "Attributes.h"
#include "FaceInfo.h"
#include "FEProblemBase.h"
#include "MooseLinearVariableFV.h"
#include "MooseMesh.h"
#include "SystemBase.h"

#include "libmesh/petsc_vector.h"

#include <algorithm>

namespace
{
/// Source of the face value of a (cell, side): no face, the interpolation of the two cell values,
/// or else a boundary entry index
constexpr int no_face = -2;
constexpr int internal_face = -1;

/// Whether the device can address host memory, as on the host backends
constexpr bool device_addresses_host =
    ::Kokkos::SpaceAccessibility<ExecSpace, ::Kokkos::HostSpace>::accessible;

/// The owned entries of the gradient component vectors, for the duration of a computation: in
/// place when the device can address them, or else a device copy written back by restore()
class GradientComponents
{
public:
  GradientComponents(FVGradientMethod::GradientContainer & gradient,
                     const dof_id_type num_local_dofs)
    : _dim(gradient.size()), _vecs(_dim), _arrays(_dim), _staged(_dim)
  {
    mooseAssert(_dim <= 3, "At most three gradient components");
    for (const auto d : make_range(_dim))
    {
      auto * const petsc_vector = dynamic_cast<libMesh::PetscVector<Number> *>(gradient[d].get());
      mooseAssert(petsc_vector, "Gradient components must be PETSc vectors");
      _vecs[d] = petsc_vector->vec();
      PetscMemType mtype;
      LibmeshPetscCallQ(VecGetArrayAndMemType(_vecs[d], &_arrays[d], &mtype));
      if (PetscMemTypeHost(mtype) && !device_addresses_host)
      {
        _staged[d].createDevice(num_local_dofs);
        _staged[d].aliasHost(_arrays[d]);
        _staged[d].copyToDevice();
        _components[d] = _staged[d].deviceData();
      }
      else
        _components[d] = _arrays[d];
    }
  }

  /// Device pointer to a component, or nullptr past the number of components
  Real * operator[](const unsigned int d) const { return _components[d]; }

  /// Write back and return the arrays; the device work on them must have finished
  void restore()
  {
    for (const auto d : make_range(_dim))
    {
      if (_staged[d].isAlloc())
        _staged[d].copyToHost();
      LibmeshPetscCallQ(VecRestoreArrayAndMemType(_vecs[d], &_arrays[d]));
    }
  }

private:
  const unsigned int _dim;
  std::vector<Vec> _vecs;
  std::vector<PetscScalar *> _arrays;
  std::vector<Moose::Kokkos::Array<Real>> _staged;
  Real * _components[3] = {nullptr, nullptr, nullptr};
};

/// The face value of a boundary side: the value relation of a boundary condition, or the cell
/// value if bc is negative
struct BoundaryFaceValue
{
  /// Compact index of the boundary condition in VariableData::bc_data
  int bc;
  /// Index of the face in the boundary condition's relation arrays
  MOOSE_KOKKOS_INDEX_TYPE index;
};
}

namespace Moose::Kokkos
{

struct GreenGaussGradient::VariableData
{
  /// Owned cells of the variable's blocks, in contiguous ID order
  Array<ContiguousElementID> cells;
  /// Whether each cell of cells is in an axisymmetric block
  Array<bool> rz;
  /// Smallest side length (Elem::hmin()) of each cell of cells, for the gradient limiter
  Array<Real> hmin;
  /// Face value source of each (cell of cells, side), with the cell index fastest
  Array2D<int> sources;
  /// Boundary face values referenced by sources
  Array<BoundaryFaceValue> boundary;
  /// Value relation arrays of the Kokkos boundary conditions of the variable
  Array<LinearFVBoundaryCondition::BoundaryData> bc_data;
  /// The Kokkos boundary conditions of the variable
  std::vector<LinearFVBoundaryCondition *> bcs;
  /// Whether any cell is in an axisymmetric block
  bool has_rz = false;
};

GreenGaussGradient::GreenGaussGradient() = default;

GreenGaussGradient::~GreenGaussGradient() = default;

GreenGaussGradient::VariableData &
GreenGaussGradient::variableData(SystemBase & system, const unsigned int variable_number)
{
  auto & data_ptr = _variable_data[{system.number(), variable_number}];
  if (data_ptr)
    return *data_ptr;
  data_ptr = std::make_unique<VariableData>();
  auto & data = *data_ptr;

  auto & problem = system.feProblem();
  auto & mesh = problem.mesh();
  const auto & kokkos_mesh = *mesh.getKokkosMesh();
  const auto & fv_mesh = problem.getKokkosFVMesh();
  const auto & variable = system.getVariable(0, variable_number);
  if (!dynamic_cast<const MooseLinearVariableFV<Real> *>(&variable))
    mooseError("Kokkos Green-Gauss gradients are computed for linear finite volume variables only; "
               "variable '",
               variable.name(),
               "' is not one.");

  problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(system.number())
      .template condition<AttribSystem>("LinearFVBoundaryCondition")
      .template condition<AttribKokkos>(true)
      .template condition<AttribThread>(0)
      .queryInto(data.bcs);
  data.bcs.erase(std::remove_if(data.bcs.begin(),
                                data.bcs.end(),
                                [variable_number](const auto * const bc)
                                { return bc->variable().var() != variable_number; }),
                 data.bcs.end());
  for (auto * const bc : data.bcs)
  {
    bc->setupFaces();
    if (!bc->hasBoundaryValue())
      mooseError("Kokkos boundary condition '",
                 bc->name(),
                 "' does not provide boundary values, which the Green-Gauss gradient of variable '",
                 variable.name(),
                 "' needs.");
  }

  // The owned cells of the variable, in contiguous ID order
  std::vector<std::pair<ContiguousElementID, const Elem *>> owned_cells;
  for (const auto * const elem : mesh.getMesh().active_local_element_ptr_range())
    if (variable.activeOnSubdomain(elem->subdomain_id()))
      owned_cells.emplace_back(kokkos_mesh.getContiguousElementID(elem), elem);
  std::sort(owned_cells.begin(), owned_cells.end());

  const auto & ranges = fv_mesh.getFaceRanges();
  auto range_of = [&ranges](const ContiguousFaceID face)
  {
    const auto it = std::upper_bound(ranges.begin(),
                                     ranges.end(),
                                     face,
                                     [](const ContiguousFaceID f, const FaceRange & r)
                                     { return f < r.begin; });
    return std::distance(ranges.begin(), it) - 1;
  };

  const auto num_cells = libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(owned_cells.size());
  const auto max_sides = fv_mesh.getMaxSidesPerElem();
  data.cells.create(num_cells);
  data.rz.create(num_cells);
  data.hmin.create(num_cells);
  data.sources.create(num_cells, max_sides);
  std::vector<BoundaryFaceValue> boundary;

  for (const auto i : make_range(num_cells))
  {
    const auto [elem, elem_ptr] = owned_cells[i];
    data.cells[i] = elem;

    const auto coord_type = mesh.getCoordSystem(elem_ptr->subdomain_id());
    if (coord_type == Moose::CoordinateSystemType::COORD_RSPHERICAL)
      mooseError("Kokkos Green-Gauss gradients do not support spherical coordinates, as "
                 "FVGreenGaussGradient. Variable '",
                 variable.name(),
                 "' needs gradients; on orthogonal meshes, set use_nonorthogonal_correction = "
                 "false on its diffusion kernels.");
    data.rz[i] = coord_type == Moose::CoordinateSystemType::COORD_RZ;
    data.hmin[i] = elem_ptr->hmin();
    data.has_rz |= data.rz[i];

    for (const auto side : make_range(max_sides))
    {
      const auto face = side < elem_ptr->n_sides() ? fv_mesh.getElemFace(elem, side)
                                                   : libMesh::DofObject::invalid_id;
      if (face == libMesh::DofObject::invalid_id)
      {
        data.sources(i, side) = no_face;
        continue;
      }

      const auto & fi = fv_mesh.getFaceInfo(face);
      const bool is_elem = fv_mesh.getElem(face) == elem;
      const auto * const other = is_elem ? fi.neighborPtr() : fi.elemPtr();
      if (other && variable.activeOnSubdomain(other->subdomain_id()))
      {
        data.sources(i, side) = internal_face;
        continue;
      }

      // A boundary of the variable: the value of the boundary condition covering the face, or
      // the cell value without one, as FVGreenGaussGradient
      const auto range = range_of(face);
      BoundaryFaceValue value{-1, 0};
      for (const auto bc_index : index_range(data.bcs))
        if (const auto offset = data.bcs[bc_index]->faceRangeOffset(range))
        {
          if (value.bc >= 0)
            mooseError("Variable '",
                       variable.name(),
                       "' has several Kokkos boundary conditions on the same face.");
          value = {libMesh::cast_int<int>(bc_index),
                   libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(*offset + face - ranges[range].begin)};
        }
      data.sources(i, side) = boundary.size();
      boundary.push_back(value);
    }
  }

  data.cells.copyToDevice();
  data.rz.copyToDevice();
  data.hmin.copyToDevice();
  data.sources.copyToDevice();
  data.boundary = boundary;
  data.bc_data.create(libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(data.bcs.size()));
  for (const auto i : index_range(data.bcs))
    data.bc_data[i] = data.bcs[i]->boundaryData();
  data.bc_data.copyToDevice();

  return data;
}

void
GreenGaussGradient::compute(SystemBase & system,
                            FVGradientMethod::GradientContainer & gradient,
                            const std::unordered_set<unsigned int> & variable_numbers)
{
  auto & problem = system.feProblem();
  const auto & fv_mesh = problem.getKokkosFVMesh();
  fv_mesh.requestSideData();

  auto & kokkos_systems = problem.getKokkosSystems();
  auto & kokkos_system = problem.getKokkosSystem(system.number());
  const auto solution_tag = problem.getVectorTagID(Moose::SOLUTION_TAG);

  kokkos_system.setActiveSolutionTags({solution_tag});
  kokkos_system.sync(MemcpyType::HOST_TO_DEVICE);
  kokkos_systems.copyToDevice();

  GradientComponents components(gradient, kokkos_system.getNumLocalDofs());

  const auto max_sides = fv_mesh.getMaxSidesPerElem();
  const auto radial = problem.mesh().getAxisymmetricRadialCoord();
  const auto sys_num = system.number();

  for (const auto variable_number : variable_numbers)
  {
    const auto & data = variableData(system, variable_number);
    // Boundary values may read the current gradients, as in FVGreenGaussGradient; the new ones are
    // written to other vectors
    for (auto * const bc : data.bcs)
    {
      bc->acquireGradients();
      bc->computeBoundaryValueData();
    }

    const auto & cells = data.cells;
    const auto & rz = data.rz;
    const auto & sources = data.sources;
    const auto & boundary = data.boundary;
    const auto & bc_data = data.bc_data;
    const bool has_rz = data.has_rz;
    const auto & systems = kokkos_systems;
    const auto var = variable_number;
    const auto mesh = fv_mesh;
    Real * const g0 = components[0];
    Real * const g1 = components[1];
    Real * const g2 = components[2];

    ::Kokkos::parallel_for(
        ::Kokkos::RangePolicy<ExecSpace>(0, cells.size()),
        KOKKOS_LAMBDA(const MOOSE_KOKKOS_INDEX_TYPE i) {
          const auto elem = cells[i];
          const auto & sys = systems[sys_num];
          const auto value = [&](const ContiguousElementID e)
          { return sys.getVectorDofValue(sys.getElemLocalDofIndex(e, 0, var), solution_tag); };
          const auto phi = value(elem);

          // Sum the faces in the local side order, so that the result does not depend on thread
          // timing
          Real3 grad;
          for (unsigned int side = 0; side < max_sides; ++side)
          {
            const auto source = sources(i, side);
            if (source == no_face)
              continue;

            Real face_value;
            if (source == internal_face)
            {
              const auto weight = mesh.getSideWeight(elem, side);
              face_value = weight * phi + (1 - weight) * value(mesh.getSideNeighbor(elem, side));
            }
            else
            {
              const auto & b = boundary[source];
              face_value = b.bc < 0 ? phi
                                    : bc_data[b.bc].value.coefficient[b.index] * phi +
                                          bc_data[b.bc].value.source[b.index];
            }
            grad += face_value * mesh.getSideAreaNormal(elem, side);
          }

          grad *= 1 / mesh.getElemVolume(elem);
          if (has_rz && rz[i])
            grad(radial) -= phi / mesh.getElemCentroid(elem)(radial);

          const auto dof = sys.getElemLocalDofIndex(elem, 0, var);
          g0[dof] = grad(0);
          if (g1)
            g1[dof] = grad(1);
          if (g2)
            g2[dof] = grad(2);
        });
  }

  ::Kokkos::fence();

  for (const auto variable_number : variable_numbers)
    for (auto * const bc : variableData(system, variable_number).bcs)
      bc->releaseGradients();

  components.restore();

  kokkos_system.sync(MemcpyType::DEVICE_TO_HOST);
  kokkos_system.clearActiveVariables();
  kokkos_system.clearActiveSolutionTags();
}

void
GreenGaussGradient::limit(SystemBase & system,
                          FVGradientMethod::GradientContainer & gradient,
                          const std::unordered_set<unsigned int> & variable_numbers,
                          const Moose::FV::GradientLimiterType limiter_type)
{
  if (limiter_type != Moose::FV::GradientLimiterType::Venkatakrishnan)
    mooseError("Kokkos Green-Gauss gradients support only the Venkatakrishnan limiter.");

  auto & problem = system.feProblem();
  const auto & fv_mesh = problem.getKokkosFVMesh();
  auto & kokkos_systems = problem.getKokkosSystems();
  auto & kokkos_system = problem.getKokkosSystem(system.number());
  const auto solution_tag = problem.getVectorTagID(Moose::SOLUTION_TAG);

  kokkos_system.setActiveSolutionTags({solution_tag});
  kokkos_system.sync(MemcpyType::HOST_TO_DEVICE);
  kokkos_systems.copyToDevice();

  GradientComponents components(gradient, kokkos_system.getNumLocalDofs());

  const auto max_sides = fv_mesh.getMaxSidesPerElem();
  const auto sys_num = system.number();
  // The host limiter limits the components up to the mesh dimension
  const unsigned int dim = problem.mesh().dimension();
  mooseAssert(dim <= gradient.size(), "Fewer gradient components than the mesh dimension");

  for (const auto variable_number : variable_numbers)
  {
    const auto & data = variableData(system, variable_number);
    const auto & cells = data.cells;
    const auto & sources = data.sources;
    const auto & hmin = data.hmin;
    const auto & systems = kokkos_systems;
    const auto var = variable_number;
    const auto mesh = fv_mesh;
    Real * const g[3] = {components[0], components[1], components[2]};

    // As ComputeLinearFVLimitedGradientThread, operation for operation: the gradient of each cell
    // is scaled by the smallest Venkatakrishnan factor over its faces, which compares the gradient
    // extrapolation to the face with the range of the values of the cell and its neighbors
    ::Kokkos::parallel_for(
        ::Kokkos::RangePolicy<ExecSpace>(0, cells.size()),
        KOKKOS_LAMBDA(const MOOSE_KOKKOS_INDEX_TYPE i) {
          const auto elem = cells[i];
          const auto & sys = systems[sys_num];
          const auto value = [&](const ContiguousElementID e)
          { return sys.getVectorDofValue(sys.getElemLocalDofIndex(e, 0, var), solution_tag); };
          const auto phi_elem = value(elem);

          Real max_value = phi_elem;
          Real min_value = phi_elem;
          for (unsigned int side = 0; side < max_sides; ++side)
            if (sources(i, side) == internal_face)
            {
              const auto phi_neighbor = value(mesh.getSideNeighbor(elem, side));
              max_value = phi_neighbor > max_value ? phi_neighbor : max_value;
              min_value = phi_neighbor < min_value ? phi_neighbor : min_value;
            }

          const auto dof = sys.getElemLocalDofIndex(elem, 0, var);
          Real3 raw_grad;
          for (unsigned int d = 0; d < dim; ++d)
            raw_grad(d) = g[d][dof];

          // A constant stencil keeps the raw gradient; the threshold is the host limiter's
          const auto range = max_value - min_value;
          if ((range < 0 ? -range : range) < 1e-14)
            return;

          const auto elem_centroid = mesh.getElemCentroid(elem);
          const Real h = hmin[i];
          const Real grad_mag = raw_grad.norm();
          // The host limiter's smoothing: 0.1 (|grad| h)^2, plus 1e-20 against a zero division
          const Real eps = 0.1 * (grad_mag * h) * (grad_mag * h) + 1e-20;
          const auto above = max_value - phi_elem;
          const auto below = min_value - phi_elem;
          const Real delta_max = (above < 0 ? -above : above) + eps;
          const Real delta_min = (below < 0 ? -below : below) + eps;

          Real alpha = 1.0;
          for (unsigned int side = 0; side < max_sides; ++side)
          {
            if (sources(i, side) != internal_face)
              continue;
            const Real delta_face =
                raw_grad * (mesh.getCentroid(mesh.getElemFace(elem, side)) - elem_centroid);
            const Real abs_delta_face = delta_face < 0 ? -delta_face : delta_face;
            const Real rf = delta_face >= 0.0 ? abs_delta_face / delta_max
                                              : abs_delta_face / delta_min;
            const Real beta = (2.0 * rf + 1.0) / (rf * (2.0 * rf + 1.0) + 1.0);
            alpha = beta < alpha ? beta : alpha;
          }

          for (unsigned int d = 0; d < dim; ++d)
            g[d][dof] = alpha * raw_grad(d);
        });
  }

  ::Kokkos::fence();
  components.restore();

  kokkos_system.clearActiveVariables();
  kokkos_system.clearActiveSolutionTags();
}

} // namespace Moose::Kokkos
