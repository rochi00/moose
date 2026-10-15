//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVKernel.h"
#include "KokkosLinearFVBoundaryCondition.h"
#include "KokkosFVDatum.h"

namespace Moose::Kokkos
{

/**
 * Base class for Kokkos linear finite volume kernels that contribute on faces (flux terms). One
 * thread handles one face and adds its contributions to the rows of both cells, so each face is
 * evaluated once. On boundary faces the contribution comes from the boundary face data populated
 * by the boundary conditions acting on this kernel's variable.
 *
 * Derived classes define computeInternalFaceContribution() and
 * computeBoundaryFaceContribution().
 *
 * With Problem/kokkos_linear_fv_assembly = gather, the face loops store the contributions per face
 * instead of adding them with atomic operations, and a loop over the cells sums every row in the
 * local side order. Each matrix and vector entry then receives one addition per loop, so the
 * assembled system is identical in every run.
 */
class LinearFVFluxKernel : public LinearFVKernel
{
public:
  static InputParameters validParams();

  LinearFVFluxKernel(const InputParameters & parameters);
  LinearFVFluxKernel(const LinearFVFluxKernel & object);

  virtual void initialSetup() override;

  /// Whether this kernel needs boundary value data from its boundary conditions
  virtual bool needsBoundaryValueData() const { return false; }
  /// Whether this kernel needs boundary normal gradient data from its boundary conditions
  virtual bool needsBoundaryNormalGradientData() const { return false; }

  /**
   * Contributions of an internal face, with the convention of the host linear FV kernels: elem and
   * neighbor are the diagonal and off-diagonal entries of the row of the elem cell, and rhs_elem
   * and rhs_neighbor the right-hand side entries of the two rows. The row of the neighbor cell is
   * the negation of the elem row with the entries swapped, which holds for fluxes that are
   * conserved across the face.
   */
  struct InternalFaceContribution
  {
    Real elem = 0;
    Real neighbor = 0;
    Real rhs_elem = 0;
    Real rhs_neighbor = 0;

    /// Diagonal entry of the row of the neighbor cell
    KOKKOS_FUNCTION Real neighborRowDiag() const { return -neighbor; }
    /// Off-diagonal entry of the row of the neighbor cell
    KOKKOS_FUNCTION Real neighborRowOff() const { return -elem; }
  };

  /// Contributions of a boundary face to the diagonal and the right-hand side of its cell's row
  struct BoundaryFaceContribution
  {
    Real diag = 0;
    Real rhs = 0;
  };

  /// Tag dispatch type for the loop over the internal faces
  struct InternalFaceLoop
  {
  };
  /// Tag dispatch type for the loop over the boundary faces
  struct BoundaryFaceLoop
  {
  };
  /// Tag dispatch type for the loop over the cells that sums the stored face contributions
  struct GatherLoop
  {
  };

  template <typename Derived>
  KOKKOS_FUNCTION void operator()(InternalFaceLoop, const ThreadID tid, const Derived & kernel) const;

  template <typename Derived>
  KOKKOS_FUNCTION void operator()(BoundaryFaceLoop, const ThreadID tid, const Derived & kernel) const;

  template <typename Derived>
  KOKKOS_FUNCTION void operator()(GatherLoop, const ThreadID tid, const Derived & kernel) const;

  /**
   * Default methods to prevent compile errors when the contributions are not defined in the
   * derived class
   */
  ///@{
  template <typename Derived>
  KOKKOS_FUNCTION InternalFaceContribution
  computeInternalFaceContribution(const FaceDatum & datum) const;
  template <typename Derived>
  KOKKOS_FUNCTION BoundaryFaceContribution computeBoundaryFaceContribution(
      const FaceDatum & datum, const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;
  ///@}

protected:
  virtual void dispatchContributions() override;

  /// Matrix positions of the contributions of an internal face of this kernel
  struct InternalFace
  {
    /// Matrix nonzero positions of the diagonal and off-diagonal entries of the elem cell row and
    /// of the neighbor cell row; -1 for a row this kernel does not assemble, because its cell is
    /// not owned or not in the kernel's blocks
    PetscInt elem_diag_slot;
    PetscInt elem_off_slot;
    PetscInt neighbor_diag_slot;
    PetscInt neighbor_off_slot;
  };

  /// A range of boundary faces of this kernel and the boundary condition that applies on it
  struct BoundaryRange
  {
    /// Whether the variable lives on the neighbor side of the faces' FaceInfo
    bool flip;
    /// Compact boundary condition index into _bc_data
    int bc_index;
    /// Index in the boundary condition's relation arrays of the first face of the range
    MOOSE_KOKKOS_INDEX_TYPE bc_offset;
  };

  /// Boundary value coefficient for this boundary face
  KOKKOS_FUNCTION Real boundaryValueCoefficient(const int bc_index,
                                                const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;
  /// Boundary value source for this boundary face
  KOKKOS_FUNCTION Real boundaryValueSource(const int bc_index,
                                           const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;
  /// Boundary normal gradient coefficient for this boundary face
  KOKKOS_FUNCTION Real boundaryNormalGradientCoefficient(
      const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;
  /// Boundary normal gradient source for this boundary face
  KOKKOS_FUNCTION Real boundaryNormalGradientSource(
      const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const;
  /// Whether the normal gradient relation of a boundary condition is a diffusive flux, see
  /// LinearFVBoundaryCondition::includesMaterialPropertyMultiplier()
  KOKKOS_FUNCTION bool boundaryIncludesMaterialPropertyMultiplier(const int bc_index) const
  {
    return _bc_includes_material_property_multiplier[bc_index];
  }
  /// Whether a boundary condition needs the non-orthogonal correction of a diffusion kernel, see
  /// LinearFVBoundaryCondition::needsBoundaryNonorthogonalCorrection()
  KOKKOS_FUNCTION bool boundaryNeedsNonorthogonalCorrection(const int bc_index) const
  {
    return _bc_needs_nonorthogonal_correction[bc_index];
  }

  /// Boundary condition data handles indexed by compact boundary condition index
  Array<LinearFVBoundaryCondition::BoundaryData> _bc_data;
  /// LinearFVBoundaryCondition::includesMaterialPropertyMultiplier() by compact boundary condition
  /// index
  Array<bool> _bc_includes_material_property_multiplier;
  /// LinearFVBoundaryCondition::needsBoundaryNonorthogonalCorrection() by compact boundary
  /// condition index
  Array<bool> _bc_needs_nonorthogonal_correction;
  /// Face ranges with the variable on both sides and at least one row this kernel assembles
  FaceRangeList _internal_ranges;
  /// Matrix positions of the contributions of each internal face, indexed by thread
  Array<InternalFace> _internal_faces;
  /// Face ranges where the variable lives on one side only and a Kokkos boundary condition applies
  FaceRangeList _boundary_ranges;
  /// Boundary condition of each range of _boundary_ranges
  Array<BoundaryRange> _boundary_range_bcs;
  /// Whether the variable lives on the neighbor side of each range of _boundary_ranges (host only)
  std::vector<bool> _host_boundary_range_flips;
  /// Matrix nonzero position of the diagonal of the cell row of each boundary face, indexed by
  /// thread
  Array<PetscInt> _boundary_diag_slots;

  /// Whether the contributions are stored per face and summed per cell (gather mode)
  bool _gather = false;
  /// Contributions stored per face in gather mode: for internal face i, the diagonal, off-diagonal
  /// and right-hand side entries of the elem cell row at 6 * i and of the neighbor cell row at
  /// 6 * i + 3; for boundary face j, the diagonal and right-hand side entries at
  /// _boundary_scratch_offset + 2 * j. Device only.
  Array<Real> _face_scratch;
  /// Start of the boundary face entries in _face_scratch
  dof_id_type _boundary_scratch_offset = 0;
  /// Cells whose rows this kernel assembles, in contiguous ID order (gather mode)
  Array<ContiguousElementID> _gather_cells;
  /// Matrix nonzero position of the diagonal of each gathered row
  Array<PetscInt> _gather_diag_slots;
  /// Position in _face_scratch of the contribution of each (gathered row, side), or -1 if the side
  /// contributes nothing; indexed with the row fastest
  Array2D<std::int64_t> _gather_sources;
  /// Matrix nonzero position of the off-diagonal entry of each (gathered row, internal side)
  Array2D<PetscInt> _gather_off_slots;
  /// Maximum number of sides of a cell
  unsigned int _max_sides = 0;

  std::unique_ptr<DispatcherBase> _internal_face_dispatcher;
  std::unique_ptr<DispatcherBase> _boundary_face_dispatcher;
  std::unique_ptr<DispatcherBase> _gather_dispatcher;
};

KOKKOS_FUNCTION inline Real
LinearFVFluxKernel::boundaryValueCoefficient(const int bc_index,
                                             const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  KOKKOS_ASSERT(bc_index >= 0);
  return _bc_data[bc_index].value.coefficient[bc_face_index];
}

KOKKOS_FUNCTION inline Real
LinearFVFluxKernel::boundaryValueSource(const int bc_index,
                                        const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  KOKKOS_ASSERT(bc_index >= 0);
  return _bc_data[bc_index].value.source[bc_face_index];
}

KOKKOS_FUNCTION inline Real
LinearFVFluxKernel::boundaryNormalGradientCoefficient(
    const int bc_index, const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  KOKKOS_ASSERT(bc_index >= 0);
  return _bc_data[bc_index].normal_gradient.coefficient[bc_face_index];
}

KOKKOS_FUNCTION inline Real
LinearFVFluxKernel::boundaryNormalGradientSource(const int bc_index,
                                                 const MOOSE_KOKKOS_INDEX_TYPE bc_face_index) const
{
  KOKKOS_ASSERT(bc_index >= 0);
  return _bc_data[bc_index].normal_gradient.source[bc_face_index];
}

template <typename Derived>
KOKKOS_FUNCTION LinearFVFluxKernel::InternalFaceContribution
LinearFVFluxKernel::computeInternalFaceContribution(const FaceDatum &) const
{
  ::Kokkos::abort("Default computeInternalFaceContribution() should never be called. Make sure "
                  "you properly redefined this method in your class without typos.");

  return {};
}

template <typename Derived>
KOKKOS_FUNCTION LinearFVFluxKernel::BoundaryFaceContribution
LinearFVFluxKernel::computeBoundaryFaceContribution(const FaceDatum &,
                                                    const int,
                                                    const MOOSE_KOKKOS_INDEX_TYPE) const
{
  ::Kokkos::abort("Default computeBoundaryFaceContribution() should never be called. Make sure "
                  "you properly redefined this method in your class without typos.");

  return {};
}

template <typename Derived>
KOKKOS_FUNCTION void
LinearFVFluxKernel::operator()(InternalFaceLoop, const ThreadID tid, const Derived & kernel) const
{
  const auto & face = _internal_faces[tid];
  const FaceDatum datum(
      _internal_ranges.face(tid, _internal_ranges.range(tid)), false, kokkosFVMesh());
  const auto contribution = kernel.template computeInternalFaceContribution<Derived>(datum);

  if (_gather)
  {
    const auto i = 6 * static_cast<std::int64_t>(tid);
    _face_scratch[i] = contribution.elem;
    _face_scratch[i + 1] = contribution.neighbor;
    _face_scratch[i + 2] = contribution.rhs_elem;
    _face_scratch[i + 3] = contribution.neighborRowDiag();
    _face_scratch[i + 4] = contribution.neighborRowOff();
    _face_scratch[i + 5] = contribution.rhs_neighbor;
    return;
  }

  KOKKOS_ASSERT(_var.components() == 1);
  const auto & sys = kokkosSystem(_var.sys());

  if (face.elem_diag_slot >= 0)
  {
    if (_assemble_matrix)
    {
      kernel.accumulateTaggedMatrixAtSlot(contribution.elem, face.elem_diag_slot);
      kernel.accumulateTaggedMatrixAtSlot(contribution.neighbor, face.elem_off_slot);
    }
    if (_assemble_rhs)
      kernel.accumulateTaggedVector(contribution.rhs_elem,
                                    sys.getElemLocalDofIndex(datum.elem(), 0, _var.var()));
  }

  if (face.neighbor_diag_slot >= 0)
  {
    if (_assemble_matrix)
    {
      kernel.accumulateTaggedMatrixAtSlot(contribution.neighborRowDiag(), face.neighbor_diag_slot);
      kernel.accumulateTaggedMatrixAtSlot(contribution.neighborRowOff(), face.neighbor_off_slot);
    }
    if (_assemble_rhs)
      kernel.accumulateTaggedVector(contribution.rhs_neighbor,
                                    sys.getElemLocalDofIndex(datum.neighbor(), 0, _var.var()));
  }
}

template <typename Derived>
KOKKOS_FUNCTION void
LinearFVFluxKernel::operator()(BoundaryFaceLoop, const ThreadID tid, const Derived & kernel) const
{
  const auto range = _boundary_ranges.range(tid);
  const auto & bc = _boundary_range_bcs[range];
  const FaceDatum datum(_boundary_ranges.face(tid, range), bc.flip, kokkosFVMesh());
  const auto contribution = kernel.template computeBoundaryFaceContribution<Derived>(
      datum, bc.bc_index, bc.bc_offset + _boundary_ranges.indexInRange(tid, range));

  if (_gather)
  {
    const auto j = _boundary_scratch_offset + 2 * static_cast<std::int64_t>(tid);
    _face_scratch[j] = contribution.diag;
    _face_scratch[j + 1] = contribution.rhs;
    return;
  }

  KOKKOS_ASSERT(_var.components() == 1);
  if (_assemble_matrix)
    kernel.accumulateTaggedMatrixAtSlot(contribution.diag, _boundary_diag_slots[tid]);
  if (_assemble_rhs)
    kernel.accumulateTaggedVector(
        contribution.rhs,
        kokkosSystem(_var.sys()).getElemLocalDofIndex(datum.elem(), 0, _var.var()));
}

template <typename Derived>
KOKKOS_FUNCTION void
LinearFVFluxKernel::operator()(GatherLoop, const ThreadID tid, const Derived & kernel) const
{
  Real diag = 0;
  Real rhs = 0;

  // Sum the row in the local side order, so that the result does not depend on thread timing
  for (unsigned int side = 0; side < _max_sides; ++side)
  {
    const auto source = _gather_sources(tid, side);
    if (source < 0)
      continue;

    diag += _face_scratch[source];
    if (source >= static_cast<std::int64_t>(_boundary_scratch_offset))
      rhs += _face_scratch[source + 1];
    else
    {
      // The off-diagonal entry of a row and side belongs to this face only
      if (_assemble_matrix)
        kernel.template accumulateTaggedMatrixAtSlot<AccumulationMode::NonAtomic>(
            _face_scratch[source + 1], _gather_off_slots(tid, side));
      rhs += _face_scratch[source + 2];
    }
  }

  KOKKOS_ASSERT(_var.components() == 1);
  if (_assemble_matrix)
    kernel.template accumulateTaggedMatrixAtSlot<AccumulationMode::NonAtomic>(
        diag, _gather_diag_slots[tid]);
  if (_assemble_rhs)
    kernel.template accumulateTaggedVector<AccumulationMode::NonAtomic>(
        rhs, kokkosSystem(_var.sys()).getElemLocalDofIndex(_gather_cells[tid], 0, _var.var()));
}

} // namespace Moose::Kokkos
