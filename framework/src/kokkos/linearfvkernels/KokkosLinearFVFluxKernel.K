//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVFluxKernel.h"
#include "KokkosLinearFVBoundaryCondition.h"
#include "LinearFVBoundaryCondition.h"

#include "Attributes.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"

#include <optional>

namespace Moose::Kokkos
{

InputParameters
LinearFVFluxKernel::validParams()
{
  auto params = LinearFVKernel::validParams();
  params.registerSystemAttributeName("LinearFVFluxKernel");
  params.addPrivateParam<bool>("_kokkos_flux_kernel", true);
  return params;
}

LinearFVFluxKernel::LinearFVFluxKernel(const InputParameters & parameters)
  : LinearFVKernel(parameters)
{
}

LinearFVFluxKernel::LinearFVFluxKernel(const LinearFVFluxKernel & object)
  : LinearFVKernel(object),
    _bc_data(object._bc_data),
    _bc_includes_material_property_multiplier(object._bc_includes_material_property_multiplier),
    _bc_needs_nonorthogonal_correction(object._bc_needs_nonorthogonal_correction),
    _internal_ranges(object._internal_ranges),
    _internal_faces(object._internal_faces),
    _boundary_ranges(object._boundary_ranges),
    _boundary_range_bcs(object._boundary_range_bcs),
    _boundary_diag_slots(object._boundary_diag_slots),
    _gather(object._gather),
    _face_scratch(object._face_scratch),
    _boundary_scratch_offset(object._boundary_scratch_offset),
    _gather_cells(object._gather_cells),
    _gather_diag_slots(object._gather_diag_slots),
    _gather_sources(object._gather_sources),
    _gather_off_slots(object._gather_off_slots),
    _max_sides(object._max_sides)
{
}

void
LinearFVFluxKernel::initialSetup()
{
  const auto sys_num = _var.sys();
  const auto var_num = _var.var();
  const auto var_name = _var.mooseVar()->name();

  // The Kokkos linear FV BCs of this kernel's variable, by compact BC index
  std::vector<LinearFVBoundaryCondition *> bcs;
  {
    std::vector<LinearFVBoundaryCondition *> all_bcs;
    _fe_problem.theWarehouse()
        .query()
        .template condition<AttribSysNum>(sys_num)
        .template condition<AttribSystem>("LinearFVBoundaryCondition")
        .template condition<AttribKokkos>(true)
        .template condition<AttribThread>(0)
        .queryInto(all_bcs);
    for (auto * bc : all_bcs)
      if (bc->variable().var() == var_num)
      {
        bc->setupFaces();
        bcs.push_back(bc);
      }
  }

  // Precompute the matrix nonzero position of every contribution so the assembly loops do not
  // search the sparsity pattern for each entry
  const auto & sys = kokkosSystem(sys_num);
  auto matrix_slot = [&](const ContiguousElementID row_elem, const ContiguousElementID col_elem)
  {
    const auto slot = sys.getMatrixSlot(sys.getElemLocalDofIndex(row_elem, 0, var_num),
                                        sys.getElemGlobalDofIndex(col_elem, 0, var_num));
    if (slot < 0)
      mooseError("Kokkos linear FV flux kernel '",
                 name(),
                 "': a matrix entry of variable '",
                 var_name,
                 "' is missing from the sparsity pattern.");
    return slot;
  };

  const auto & fv_mesh = kokkosFVMesh();
  const auto num_local_elems = fv_mesh.getNumLocalElements();
  const auto * const variable = _var.mooseFieldVar();
  const auto & mesh_ranges = fv_mesh.getFaceRanges();

  // Whether the variable lives on a block
  auto variable_on = [&](const SubdomainID subdomain)
  { return subdomain != Moose::INVALID_BLOCK_ID && variable->activeOnSubdomain(subdomain); };

  std::vector<FaceRange> internal_ranges, boundary_ranges;
  std::vector<BoundaryRange> boundary_range_bcs;
  // Rows assembled on the elem and neighbor side of each internal range
  std::vector<std::pair<bool, bool>> internal_rows;

  for (const auto r : index_range(mesh_ranges))
  {
    const auto & range = mesh_ranges[r];
    const bool on_elem = variable_on(range.elem_subdomain);
    const bool on_neighbor = variable_on(range.neighbor_subdomain);
    // This kernel assembles the row of a cell that is owned and in the kernel's blocks
    const bool elem_row = range.elem_owned && hasBlocks(range.elem_subdomain);
    const bool neighbor_row = range.neighbor_owned && hasBlocks(range.neighbor_subdomain);

    if (on_elem && on_neighbor)
    {
      if (elem_row || neighbor_row)
      {
        internal_ranges.push_back(range);
        internal_rows.emplace_back(elem_row, neighbor_row);
      }
      continue;
    }

    // A boundary for the variable; it contributes through a boundary condition, if one applies
    const bool flip = on_neighbor;
    if (!(flip ? neighbor_row : elem_row))
      continue;

    const auto & range_boundaries = fv_mesh.getFaceRangeBoundaryIDs(r, flip);
    if (range_boundaries.size() > 1)
      mooseError("Kokkos linear FV flux kernel '",
                 name(),
                 "': boundary faces of variable '",
                 var_name,
                 "' carry several boundary IDs ",
                 Moose::stringify(range_boundaries),
                 ". Multiple boundary conditions on the same face are not supported.");

    std::optional<BoundaryRange> range_bc;
    for (const auto bc_index : index_range(bcs))
      if (const auto offset = bcs[bc_index]->faceRangeOffset(r))
      {
        // Ambiguous boundary data would make flux kernels non-deterministic.
        if (range_bc)
          mooseError("Kokkos linear FV flux kernel '",
                     name(),
                     "' found multiple Kokkos boundary conditions for variable '",
                     var_name,
                     "' on boundary ",
                     range_boundaries.front(),
                     ".");
        range_bc = BoundaryRange{flip, libMesh::cast_int<int>(bc_index), *offset};
      }

    if (range_bc)
    {
      boundary_ranges.push_back(range);
      boundary_range_bcs.push_back(*range_bc);
    }
  }

  // Matrix positions per face, in thread order
  std::vector<InternalFace> internal_faces;
  for (const auto i : index_range(internal_ranges))
    for (const auto face : make_range(internal_ranges[i].begin, internal_ranges[i].end))
    {
      const auto elem = fv_mesh.getElem(face);
      const auto neighbor = fv_mesh.getNeighbor(face);
      const auto [elem_row, neighbor_row] = internal_rows[i];
      internal_faces.push_back({elem_row ? matrix_slot(elem, elem) : -1,
                                elem_row ? matrix_slot(elem, neighbor) : -1,
                                neighbor_row ? matrix_slot(neighbor, neighbor) : -1,
                                neighbor_row ? matrix_slot(neighbor, elem) : -1});
    }

  std::vector<PetscInt> boundary_diag_slots;
  for (const auto i : index_range(boundary_ranges))
    for (const auto face : make_range(boundary_ranges[i].begin, boundary_ranges[i].end))
    {
      const auto elem =
          boundary_range_bcs[i].flip ? fv_mesh.getNeighbor(face) : fv_mesh.getElem(face);
      boundary_diag_slots.push_back(matrix_slot(elem, elem));
    }

  _internal_ranges.create(internal_ranges);
  _internal_faces = internal_faces;
  _boundary_ranges.create(boundary_ranges);
  _boundary_range_bcs = boundary_range_bcs;
  _boundary_diag_slots = boundary_diag_slots;

  _bc_data.create(libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(bcs.size()));
  _bc_includes_material_property_multiplier.create(
      libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(bcs.size()));
  _bc_needs_nonorthogonal_correction.create(libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(bcs.size()));
  for (const auto i : index_range(bcs))
  {
    _bc_data[i] = bcs[i]->boundaryData();
    _bc_includes_material_property_multiplier[i] = bcs[i]->includesMaterialPropertyMultiplier();
    _bc_needs_nonorthogonal_correction[i] = bcs[i]->needsBoundaryNonorthogonalCorrection();
  }
  _bc_data.copyToDevice();
  _bc_includes_material_property_multiplier.copyToDevice();
  _bc_needs_nonorthogonal_correction.copyToDevice();

  _gather = _fe_problem.getParam<MooseEnum>("kokkos_linear_fv_assembly") == "gather";
  if (_gather)
  {
    _max_sides = fv_mesh.getMaxSidesPerElem();
    _boundary_scratch_offset = 6 * internal_faces.size();
    _face_scratch.createDevice(_boundary_scratch_offset + 2 * boundary_diag_slots.size());

    // Scratch position and off-diagonal matrix position of every (owned cell, side) this kernel
    // assembles
    const auto num_entries = std::size_t(num_local_elems) * _max_sides;
    std::vector<std::int64_t> sources(num_entries, -1);
    std::vector<PetscInt> off_slots(num_entries, -1);
    auto add_source = [&](const ContiguousElementID elem,
                          const unsigned int side,
                          const std::int64_t source,
                          const PetscInt off_slot)
    {
      const auto entry = std::size_t(elem) * _max_sides + side;
      mooseAssert(sources[entry] < 0, "A cell side must have one face");
      sources[entry] = source;
      off_slots[entry] = off_slot;
    };

    std::int64_t i = 0;
    for (const auto & range : internal_ranges)
      for (const auto face : make_range(range.begin, range.end))
      {
        const auto & slots = internal_faces[i];
        if (slots.elem_diag_slot >= 0)
          add_source(fv_mesh.getElem(face), fv_mesh.getElemSide(face), 6 * i, slots.elem_off_slot);
        if (slots.neighbor_diag_slot >= 0)
          add_source(fv_mesh.getNeighbor(face),
                     fv_mesh.getNeighborSide(face),
                     6 * i + 3,
                     slots.neighbor_off_slot);
        ++i;
      }
    std::int64_t j = 0;
    for (const auto r : index_range(boundary_ranges))
      for (const auto face : make_range(boundary_ranges[r].begin, boundary_ranges[r].end))
      {
        const bool flip = boundary_range_bcs[r].flip;
        add_source(flip ? fv_mesh.getNeighbor(face) : fv_mesh.getElem(face),
                   flip ? fv_mesh.getNeighborSide(face) : fv_mesh.getElemSide(face),
                   _boundary_scratch_offset + 2 * j++,
                   -1);
      }

    std::vector<ContiguousElementID> gather_cells;
    for (const auto elem : make_range(num_local_elems))
      for (const auto side : make_range(_max_sides))
        if (sources[std::size_t(elem) * _max_sides + side] >= 0)
        {
          gather_cells.push_back(elem);
          break;
        }

    const auto num_rows = libMesh::cast_int<MOOSE_KOKKOS_INDEX_TYPE>(gather_cells.size());
    _gather_cells = gather_cells;
    _gather_diag_slots.create(num_rows);
    _gather_sources.create(num_rows, _max_sides);
    _gather_off_slots.create(num_rows, _max_sides);
    for (const auto row : make_range(num_rows))
    {
      const auto elem = gather_cells[row];
      _gather_diag_slots[row] = matrix_slot(elem, elem);
      for (const auto side : make_range(_max_sides))
      {
        _gather_sources(row, side) = sources[std::size_t(elem) * _max_sides + side];
        _gather_off_slots(row, side) = off_slots[std::size_t(elem) * _max_sides + side];
      }
    }
    _gather_diag_slots.copyToDevice();
    _gather_sources.copyToDevice();
    _gather_off_slots.copyToDevice();
  }

  // Non-Kokkos BCs are never called by Kokkos flux kernels; they would be silently
  // ignored, producing wrong results. Reject them early for this variable.
  std::vector<::LinearFVBoundaryCondition *> non_kokkos_bcs;
  _fe_problem.theWarehouse()
      .query()
      .template condition<AttribSysNum>(sys_num)
      .template condition<AttribSystem>("LinearFVBoundaryCondition")
      .template condition<AttribKokkos>(false)
      .template condition<AttribThread>(0)
      .queryInto(non_kokkos_bcs);

  for (auto * bc : non_kokkos_bcs)
  {
    const auto & variable = bc->variable();

    if (variable.number() == var_num)
      mooseError("Kokkos linear FV flux kernel '",
                 name(),
                 "' found non-Kokkos boundary condition '",
                 bc->name(),
                 "' (type '",
                 bc->type(),
                 "') for variable '",
                 variable.name(),
                 "'. Non-Kokkos BCs are not called by Kokkos flux kernels and will produce "
                 "incorrect results. Use the Kokkos equivalent, e.g. "
                 "'KokkosLinearFVFunctorDirichletBC'.");
  }
}

void
LinearFVFluxKernel::dispatchContributions()
{
  if (_internal_ranges.size())
  {
    Policy policy(0, _internal_ranges.size());
    if (!_internal_face_dispatcher)
      _internal_face_dispatcher = DispatcherRegistry::build<InternalFaceLoop>(this, type());
    _internal_face_dispatcher->parallelForAsync(policy);
  }

  if (_boundary_ranges.size())
  {
    Policy policy(0, _boundary_ranges.size());
    if (!_boundary_face_dispatcher)
      _boundary_face_dispatcher = DispatcherRegistry::build<BoundaryFaceLoop>(this, type());
    _boundary_face_dispatcher->parallelForAsync(policy);
  }

  if (_gather && _gather_cells.size())
  {
    Policy policy(0, _gather_cells.size());
    if (!_gather_dispatcher)
      _gather_dispatcher = DispatcherRegistry::build<GatherLoop>(this, type());
    _gather_dispatcher->parallelForAsync(policy);
  }
}

} // namespace Moose::Kokkos
