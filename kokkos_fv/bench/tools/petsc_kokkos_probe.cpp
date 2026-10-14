// Answers design question Q8: does PETSc's Kokkos backend run on Kokkos's default execution space
// instance? If it does, work that MOOSE dispatches on the default instance is ordered before
// PETSc's MatSetValuesCOO/VecGhostUpdate on the same stream, and the fence before handing the
// matrix to PETSc is not needed for ordering.
//
// Build (in the MOOSE conda env, or against the GPU PETSc):
//   mpicxx -std=c++17 -O2 [-fopenmp] petsc_kokkos_probe.cpp -I$PETSC_DIR/include -L$PETSC_DIR/lib \
//     -lpetsc -lkokkoscore -Wl,-rpath,$PETSC_DIR/lib -o petsc_kokkos_probe
// Run: ./petsc_kokkos_probe   (PETSc initializes Kokkos when the first Kokkos object is created)

#include <petscvec.h>
#include <petsc_kokkos.hpp>
#include <Kokkos_Core.hpp>

#include <cstdio>

int
main(int argc, char ** argv)
{
  PetscCall(PetscInitialize(&argc, &argv, nullptr, nullptr));
  {
    Vec v;
    PetscCall(VecCreate(PETSC_COMM_WORLD, &v));
    PetscCall(VecSetSizes(v, 16, PETSC_DECIDE));
    PetscCall(VecSetType(v, VECKOKKOS));
    PetscCall(VecSet(v, 1.0));

    const auto petsc_exec = PetscGetKokkosExecutionSpace();
    const Kokkos::DefaultExecutionSpace default_exec;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                          "execution space: %s\n"
                          "PETSc instance id: %u\n"
                          "Kokkos default instance id: %u\n"
                          "same instance: %s\n",
                          Kokkos::DefaultExecutionSpace::name(),
                          petsc_exec.impl_instance_id(),
                          default_exec.impl_instance_id(),
                          petsc_exec.impl_instance_id() == default_exec.impl_instance_id() ? "yes"
                                                                                          : "no"));
    PetscCall(VecDestroy(&v));
  }
  PetscCall(PetscFinalize());
  return 0;
}
