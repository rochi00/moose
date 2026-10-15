//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosFVSystemAssembly.h"
#include "KokkosSystem.h"

#include "SystemBase.h"

#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"

namespace Moose::Kokkos
{

namespace
{
/// Whether the device can address host memory, as on the host backends
constexpr bool device_addresses_host =
    ::Kokkos::SpaceAccessibility<ExecSpace, ::Kokkos::HostSpace>::accessible;

/// Zero a device array without waiting for it to finish
void
zeroAsync(const Array<Real> & array)
{
  Real * const values = array.deviceData();
  ::Kokkos::parallel_for(
      ::Kokkos::RangePolicy<ExecSpace>(0, array.size()),
      KOKKOS_LAMBDA(const std::size_t i) { values[i] = 0; });
}
}

FVSystemAssembly::FVSystemAssembly(SystemBase & system, const System & kokkos_system)
  : _system(system),
    _kokkos_system(kokkos_system),
    _num_owned_rows(kokkos_system.getNumLocalDofs()),
    // The sparsity pattern lists the owned rows first
    _num_owned_nonzeros(kokkos_system.getSparsity().row_ptr[_num_owned_rows])
{
}

void
FVSystemAssembly::begin(const std::set<TagID> & vector_tags, const std::set<TagID> & matrix_tags)
{
  _active_matrix_tags.clear();
  for (const auto tag : matrix_tags)
  {
    if (!_system.hasMatrix(tag))
      continue;

    auto * const petsc_matrix =
        dynamic_cast<libMesh::PetscMatrix<PetscScalar> *>(&_system.getMatrix(tag));
    mooseAssert(petsc_matrix, "Kokkos linear FV matrices must be PETSc matrices");
    const auto mat = petsc_matrix->mat();

    auto & matrix = _matrices[tag];
    if (matrix.mat != mat)
    {
      PetscMemType mtype;
      LibmeshPetscCallQ(MatGetCurrentMemType(mat, &mtype));
      matrix.is_host = PetscMemTypeHost(mtype);

      // PETSc sorts and biases the index arrays in place, so they are filled for each matrix
      std::vector<PetscInt> coo_i, coo_j;
      _kokkos_system.buildCooIndices(coo_i, coo_j);
      LibmeshPetscCallQ(
          MatSetPreallocationCOO(mat, _num_owned_nonzeros, coo_i.data(), coo_j.data()));

      if (matrix.is_host && !device_addresses_host)
        matrix.values.create(_num_owned_nonzeros);
      else
        matrix.values.createDevice(_num_owned_nonzeros);
      matrix.mat = mat;
    }

    zeroAsync(matrix.values);
    _active_matrix_tags.insert(tag);
  }

  for (const auto tag : vector_tags)
  {
    if (!_system.hasVector(tag))
      continue;

    auto * const petsc_vector =
        dynamic_cast<libMesh::PetscVector<PetscScalar> *>(&_system.getVector(tag));
    mooseAssert(petsc_vector, "Kokkos linear FV vectors must be PETSc vectors");

    auto & vector = _vectors[tag];
    vector.vec = petsc_vector->vec();
    PetscMemType mtype;
    LibmeshPetscCallQ(VecGetArrayAndMemType(vector.vec, &vector.array, &mtype));

    // The objects add into the owned entries of the local array, in place when the device can
    // address it, or else into a device copy
    vector.staged = PetscMemTypeHost(mtype) && !device_addresses_host;
    if (vector.staged)
    {
      if (!vector.values.isAlloc())
        vector.values.createDevice(_num_owned_rows);
      vector.values.aliasHost(vector.array);
      vector.values.copyToDevice();
    }
    else
    {
      if (!vector.values.isAlloc())
        vector.values.init(_num_owned_rows);
      vector.values.aliasDevice(vector.array);
    }
  }
}

Real *
FVSystemAssembly::vectorValues(const TagID tag) const
{
  const auto it = _vectors.find(tag);
  return it != _vectors.end() && it->second.vec ? it->second.values.deviceData() : nullptr;
}

Real *
FVSystemAssembly::matrixValues(const TagID tag) const
{
  return _active_matrix_tags.count(tag) ? libmesh_map_find(_matrices, tag).values.deviceData()
                                        : nullptr;
}

void
FVSystemAssembly::end()
{
  for (const auto tag : _active_matrix_tags)
  {
    auto & matrix = libmesh_map_find(_matrices, tag);
    const bool to_host = matrix.is_host && !device_addresses_host;
    if (to_host)
      matrix.values.copyToHost();
    LibmeshPetscCallQ(MatSetValuesCOO(matrix.mat,
                                      to_host ? matrix.values.hostData()
                                              : matrix.values.deviceData(),
                                      ADD_VALUES));
  }

  for (auto & [tag, vector] : _vectors)
  {
    if (!vector.vec)
      continue;
    if (vector.staged)
      vector.values.copyToHost();
    LibmeshPetscCallQ(VecRestoreArrayAndMemType(vector.vec, &vector.array));
    vector.vec = nullptr;
  }

  // PETSc reads the matrix values on its own execution space instance (design question Q8); wait
  // for it before the next assembly zeroes them
  ::Kokkos::fence();
}

} // namespace Moose::Kokkos
