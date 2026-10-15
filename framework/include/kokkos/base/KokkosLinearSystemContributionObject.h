//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosDispatcher.h"
#include "KokkosMesh.h"
#include "KokkosFVMesh.h"
#include "KokkosFVSystemAssembly.h"
#include "KokkosLinearFVPropertyInterface.h"
#include "KokkosLinearFVVariableInterface.h"
#include "KokkosSystem.h"
#include "KokkosVariable.h"

#include "MooseObject.h"
#include "SetupInterface.h"
#include "FunctionInterface.h"
#include "UserObjectInterface.h"
#include "TransientInterface.h"
#include "PostprocessorInterface.h"
#include "VectorPostprocessorInterface.h"
#include "Restartable.h"
#include "MeshChangedInterface.h"
#include "TaggingInterface.h"
#include "LinearSystemContributionObject.h"

class InputParameters;
class FEProblemBase;

namespace Moose::Kokkos
{

/**
 * Base class for Kokkos objects that contribute to a linear system, i.e. the linear finite volume
 * kernels and boundary conditions. It provides the device-side helpers for accumulating into the
 * tagged vectors and matrices.
 */
class LinearSystemContributionObject : public MooseObject,
                                       public SetupInterface,
                                       public FunctionInterface,
                                       public UserObjectInterface,
                                       public TransientInterface,
                                       public PostprocessorInterface,
                                       public VectorPostprocessorInterface,
                                       public Restartable,
                                       public MeshChangedInterface,
                                       public TaggingInterface,
                                       public MeshHolder,
                                       public FVMeshHolder,
                                       public SystemHolder,
                                       public LinearFVVariableInterface,
                                       public LinearFVPropertyInterface
{
public:
  static InputParameters validParams();

  LinearSystemContributionObject(const InputParameters & parameters);
  LinearSystemContributionObject(const LinearSystemContributionObject & object);

  /**
   * Get the Kokkos variable this object contributes to
   * @returns The Kokkos variable
   */
  Variable variable() const { return _var; }

  /**
   * Compute the contributions of this object to the tagged matrices, the tagged vectors, or both,
   * in one pass. The work is dispatched without waiting for it to finish; the caller fences before
   * reading the results.
   * @param matrix Whether to add the matrix contributions
   * @param rhs Whether to add the right-hand side contributions
   * @param assembly The arrays of the tagged matrices and vectors being assembled
   */
  void computeContributions(const bool matrix, const bool rhs, const FVSystemAssembly & assembly);

protected:
  /// Whether tagged vector/matrix accumulation must use atomics
  enum class AccumulationMode
  {
    Atomic,
    NonAtomic
  };

  /// Add a contribution to a row of every tagged vector this object assembles
  template <AccumulationMode mode = AccumulationMode::Atomic>
  KOKKOS_FUNCTION void accumulateTaggedVector(Real value, dof_id_type row) const;
  /// Add a contribution to every tagged matrix this object assembles, at a precomputed nonzero
  /// position (see System::getMatrixSlot()), which avoids searching the row for the column
  template <AccumulationMode mode = AccumulationMode::Atomic>
  KOKKOS_FUNCTION void accumulateTaggedMatrixAtSlot(Real value, PetscInt slot) const;

  /// Reference to the finite element problem
  FEProblemBase & _fe_problem;

  /// Dispatch the loops of this object for computeContributions()
  virtual void dispatchContributions() = 0;

  /// Kokkos variable
  Variable _var;

  /// Whether the current computeContributions() adds the matrix contributions
  bool _assemble_matrix = false;
  /// Whether the current computeContributions() adds the right-hand side contributions
  bool _assemble_rhs = false;

  /// Most tagged matrices or vectors an object adds to in one assembly; a linear system normally
  /// assembles one of each
  static constexpr unsigned int max_tagged_targets = 4;
  /// Device arrays of the tagged vectors this object adds to in the current assembly
  Real * _vector_targets[max_tagged_targets] = {};
  /// Number of entries of _vector_targets
  unsigned int _num_vector_targets = 0;
  /// Device value arrays of the tagged matrices this object adds to in the current assembly
  Real * _matrix_targets[max_tagged_targets] = {};
  /// Number of entries of _matrix_targets
  unsigned int _num_matrix_targets = 0;

  /**
   * TODO: Move to TransientInterface
   */
  ///@{
  /// Current time
  Scalar<Real> _t;
  ///@}
};

template <LinearSystemContributionObject::AccumulationMode mode>
KOKKOS_FUNCTION inline void
LinearSystemContributionObject::accumulateTaggedVector(const Real value,
                                                       const dof_id_type row) const
{
  if (!value)
    return;

  for (unsigned int i = 0; i < _num_vector_targets; ++i)
  {
    if constexpr (mode == AccumulationMode::Atomic)
      ::Kokkos::atomic_add(&_vector_targets[i][row], value);
    else
      _vector_targets[i][row] += value;
  }
}

template <LinearSystemContributionObject::AccumulationMode mode>
KOKKOS_FUNCTION inline void
LinearSystemContributionObject::accumulateTaggedMatrixAtSlot(const Real value,
                                                             const PetscInt slot) const
{
  if (!value)
    return;

  KOKKOS_ASSERT(slot >= 0);
  for (unsigned int i = 0; i < _num_matrix_targets; ++i)
  {
    if constexpr (mode == AccumulationMode::Atomic)
      ::Kokkos::atomic_add(&_matrix_targets[i][slot], value);
    else
      _matrix_targets[i][slot] += value;
  }
}

} // namespace Moose::Kokkos
