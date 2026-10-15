//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVVariableInterface.h"

#include "FEProblemBase.h"
#include "LinearFVGradientReader.h"
#include "MooseLinearVariableFV.h"
#include "MooseObject.h"
#include "SystemBase.h"

#include "libmesh/petsc_vector.h"

namespace Moose::Kokkos
{

LinearFVVariableInterface::LinearFVVariableInterface(const MooseObject & moose_object,
                                                     FEProblemBase & problem,
                                                     const unsigned int own_system)
  : _lfv_moose_object(moose_object), _lfv_problem(problem), _read_systems({own_system})
{
}

LinearFVVariableInterface::LinearFVVariableInterface(const LinearFVVariableInterface & object)
  : _lfv_moose_object(object._lfv_moose_object), _lfv_problem(object._lfv_problem)
{
}

LinearFVCellValues
LinearFVVariableInterface::getLinearFVCellValues(const VariableName & name,
                                                 const unsigned int state)
{
  auto & variable = _lfv_problem.getActualFieldVariable(0, name);
  if (!dynamic_cast<MooseLinearVariableFV<Real> *>(&variable))
    _lfv_moose_object.mooseError(
        "Variable '", name, "' must be a linear finite volume variable (MooseLinearVariableFVReal).");

  static const std::vector<TagName> state_tags = {
      Moose::SOLUTION_TAG, Moose::OLD_SOLUTION_TAG, Moose::OLDER_SOLUTION_TAG};
  if (state >= state_tags.size())
    _lfv_moose_object.mooseError("Solution state ",
                                 state,
                                 " of variable '",
                                 name,
                                 "' is not available; the states are 0 (current), 1 (old) and 2 "
                                 "(older).");

  // Make the system keep the state, so that its tagged vector exists and is synchronized
  auto & system = variable.sys();
  system.needSolutionState(state);
  _read_systems.insert(system.number());

  return {system.number(), variable.number(), _lfv_problem.getVectorTagID(state_tags[state])};
}

void
LinearFVVariableInterface::requestLinearFVCellGradient(const VariableName & name,
                                                       LinearFVCellGradient & gradient,
                                                       const GradientMethodName & method_name)
{
  auto * const variable =
      dynamic_cast<MooseLinearVariableFV<Real> *>(&_lfv_problem.getActualFieldVariable(0, name));
  if (!variable)
    _lfv_moose_object.mooseError(
        "Variable '", name, "' must be a linear finite volume variable (MooseLinearVariableFVReal).");

  gradient.sys = variable->sys().number();
  gradient.var = variable->number();
  _read_systems.insert(gradient.sys);
  const auto & reader = method_name.empty() ? variable->requestCellGradients()
                                             : variable->requestCellGradients(method_name);
  _gradients.push_back({&gradient, &reader, {}, {}, {}, {}});
}

void
LinearFVVariableInterface::acquireGradients()
{
  // Whether the device can address host memory, as on the host backends
  constexpr bool device_addresses_host =
      ::Kokkos::SpaceAccessibility<ExecSpace, ::Kokkos::HostSpace>::accessible;

  if (_gradient_acquisitions++)
    return;

  for (auto & access : _gradients)
  {
    const auto & components = access.reader->components();
    const auto dim = components.size();
    access.global.resize(dim);
    access.local.resize(dim);
    access.arrays.resize(dim);
    access.staged.resize(dim);

    for (const auto d : make_range(dim))
    {
      auto * const petsc_vector = dynamic_cast<libMesh::PetscVector<Number> *>(components[d]);
      mooseAssert(petsc_vector, "Gradient components must be PETSc vectors");
      access.global[d] = petsc_vector->vec();
      // The local form holds the ghost entries after the owned ones, as the Kokkos system's local
      // degrees of freedom
      LibmeshPetscCallQ(VecGhostGetLocalForm(access.global[d], &access.local[d]));
      const auto vec = access.local[d] ? access.local[d] : access.global[d];

      PetscMemType mtype;
      LibmeshPetscCallQ(VecGetArrayReadAndMemType(vec, &access.arrays[d], &mtype));
      if (PetscMemTypeHost(mtype) && !device_addresses_host)
      {
        PetscInt size;
        LibmeshPetscCallQ(VecGetLocalSize(vec, &size));
        if (!access.staged[d].isAlloc())
          access.staged[d].createDevice(size);
        access.staged[d].aliasHost(const_cast<PetscScalar *>(access.arrays[d]));
        access.staged[d].copyToDevice();
        access.handle->components[d] = access.staged[d].deviceData();
      }
      else
        access.handle->components[d] = access.arrays[d];
    }
  }
}

void
LinearFVVariableInterface::releaseGradients()
{
  mooseAssert(_gradient_acquisitions, "releaseGradients() without acquireGradients()");
  if (--_gradient_acquisitions)
    return;

  for (auto & access : _gradients)
    for (const auto d : index_range(access.arrays))
    {
      const auto vec = access.local[d] ? access.local[d] : access.global[d];
      LibmeshPetscCallQ(VecRestoreArrayReadAndMemType(vec, &access.arrays[d]));
      LibmeshPetscCallQ(VecGhostRestoreLocalForm(access.global[d], &access.local[d]));
      access.handle->components[d] = nullptr;
    }
}

} // namespace Moose::Kokkos
