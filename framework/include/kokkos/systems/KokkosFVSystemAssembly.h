//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosArray.h"
#include "MooseTypes.h"

#include <petscmat.h>
#include <petscvec.h>

#include <map>
#include <set>

class SystemBase;

namespace Moose::Kokkos
{

class System;

/**
 * The tagged matrices and right-hand side vectors of one assembly of a linear system by the
 * Kokkos linear finite volume objects.
 *
 * The objects add their contributions into the arrays this class hands out: one value array per
 * matrix tag, laid out as the owned rows of the Kokkos system's sparsity pattern so that the
 * positions from System::getMatrixSlot() apply, and the local array of each tagged PETSc vector.
 * The objects only write rows of owned cells, so the matrix pattern is preallocated with the owned
 * rows only, which leaves no entries to communicate in MatSetValuesCOO(). Neither the tagged
 * vectors nor the tagged matrices go through the Kokkos system's synchronization, so an assembly
 * allocates nothing and does not wait for the device until end().
 */
class FVSystemAssembly
{
public:
  /**
   * Constructor
   * @param system The MOOSE linear system
   * @param kokkos_system The Kokkos system of the linear system
   */
  FVSystemAssembly(SystemBase & system, const System & kokkos_system);

  /**
   * Prepare an assembly: zero the matrix value arrays and get the vector arrays. The work is
   * dispatched without waiting for it to finish.
   * @param vector_tags The vector tags to assemble; tags without a vector in the system are ignored
   * @param matrix_tags The matrix tags to assemble; tags without a matrix in the system are ignored
   */
  void begin(const std::set<TagID> & vector_tags, const std::set<TagID> & matrix_tags);

  /**
   * Get the device array into which the contributions to a tagged vector are added, indexed by
   * local degree of freedom
   * @returns The array, or nullptr if the tag is not assembled
   */
  Real * vectorValues(const TagID tag) const;
  /**
   * Get the device array into which the contributions to a tagged matrix are added, indexed by
   * System::getMatrixSlot()
   * @returns The array, or nullptr if the tag is not assembled
   */
  Real * matrixValues(const TagID tag) const;

  /**
   * Hand the assembled values to PETSc. The Kokkos work that adds to the arrays must have finished.
   */
  void end();

private:
  /// The value array of a tagged matrix
  struct MatrixValues
  {
    /// The PETSc matrix the pattern was preallocated for
    Mat mat = nullptr;
    /// Whether the PETSc matrix keeps its values in host memory
    bool is_host = false;
    /// The values, on device; also on host if the device cannot address host memory and the
    /// matrix is a host matrix
    Array<Real> values;
  };

  /// The array of a tagged vector during an assembly
  struct VectorValues
  {
    /// The PETSc vector
    Vec vec = nullptr;
    /// Its local array
    PetscScalar * array = nullptr;
    /// Whether the array is in host memory that the device cannot address, so that the values
    /// are staged on device
    bool staged = false;
    /// The values the objects add into: aliases the local array, or the staging copy
    Array<Real> values;
  };

  /// The MOOSE linear system
  SystemBase & _system;
  /// The Kokkos system of the linear system
  const System & _kokkos_system;
  /// Number of owned rows
  PetscInt _num_owned_rows;
  /// Number of nonzeros of the owned rows
  PetscCount _num_owned_nonzeros;

  /// Value arrays by matrix tag, kept between assemblies
  std::map<TagID, MatrixValues> _matrices;
  /// Arrays by vector tag, for the current assembly
  std::map<TagID, VectorValues> _vectors;
  /// Matrix tags of the current assembly
  std::set<TagID> _active_matrix_tags;
};

} // namespace Moose::Kokkos
