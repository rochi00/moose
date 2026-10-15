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
#include "KokkosFVDatum.h"

namespace Moose::Kokkos
{

/**
 * Base class for Kokkos linear finite volume kernels that contribute on cells (volumetric terms).
 * Derived classes define computeCellContribution(), which returns the contributions of a cell to
 * the diagonal and the right-hand side of its own row.
 */
class LinearFVElementalKernel : public LinearFVKernel
{
public:
  static InputParameters validParams();

  LinearFVElementalKernel(const InputParameters & parameters);
  LinearFVElementalKernel(const LinearFVElementalKernel & object);

  virtual void initialSetup() override;

  /// Contributions of a cell to the diagonal and the right-hand side of its row
  struct CellContribution
  {
    Real diag = 0;
    Real rhs = 0;
  };

  /// Tag dispatch type for the loop over the cells
  struct CellLoop
  {
  };

  /**
   * Cell loop body; adds the contributions of one cell
   * @param tid The thread ID of the current cell
   * @param kernel The derived kernel object whose contribution method is invoked
   */
  template <typename Derived>
  KOKKOS_FUNCTION void operator()(CellLoop, const ThreadID tid, const Derived & kernel) const;

  /**
   * Default method to prevent compile errors when the contribution is not defined in the derived
   * class
   */
  template <typename Derived>
  KOKKOS_FUNCTION CellContribution computeCellContribution(const CellDatum &) const
  {
    ::Kokkos::abort("Default computeCellContribution() should never be called. Make sure you "
                    "properly redefined this method in your class without typos.");

    return {};
  }

protected:
  virtual void dispatchContributions() override;

  /// Matrix nonzero position of the diagonal of each block element, indexed like the dispatch
  /// thread, so the cell loop does not search the sparsity pattern
  Array<PetscInt> _diagonal_matrix_slots;

  /// Dispatcher of the cell loop
  std::unique_ptr<DispatcherBase> _cell_dispatcher;
};

template <typename Derived>
KOKKOS_FUNCTION void
LinearFVElementalKernel::operator()(CellLoop, const ThreadID tid, const Derived & kernel) const
{
  const auto elem = kokkosBlockElementID(tid);
  const CellDatum datum(elem, kokkosFVMesh());
  const auto contribution = kernel.template computeCellContribution<Derived>(datum);

  KOKKOS_ASSERT(_var.components() == 1);
  if (_assemble_matrix)
    kernel.template accumulateTaggedMatrixAtSlot<AccumulationMode::NonAtomic>(
        contribution.diag, _diagonal_matrix_slots[tid]);
  if (_assemble_rhs)
    kernel.template accumulateTaggedVector<AccumulationMode::NonAtomic>(
        contribution.rhs, kokkosSystem(_var.sys()).getElemLocalDofIndex(elem, 0, _var.var()));
}

} // namespace Moose::Kokkos
