//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVElementalKernel.h"

namespace Moose::Kokkos
{

InputParameters
LinearFVElementalKernel::validParams()
{
  auto params = LinearFVKernel::validParams();
  params.registerSystemAttributeName("LinearFVElementalKernel");
  return params;
}

LinearFVElementalKernel::LinearFVElementalKernel(const InputParameters & parameters)
  : LinearFVKernel(parameters)
{
}

LinearFVElementalKernel::LinearFVElementalKernel(const LinearFVElementalKernel & object)
  : LinearFVKernel(object), _diagonal_matrix_slots(object._diagonal_matrix_slots)
{
}

void
LinearFVElementalKernel::initialSetup()
{
  LinearFVKernel::initialSetup();

  const auto & sys = kokkosSystem(_var.sys());
  const auto var_num = _var.var();
  std::vector<PetscInt> slots(numKokkosBlockElements());
  for (const auto tid : index_range(slots))
  {
    const auto elem = kokkosBlockElementID(tid);
    slots[tid] = sys.getMatrixSlot(sys.getElemLocalDofIndex(elem, 0, var_num),
                                   sys.getElemGlobalDofIndex(elem, 0, var_num));
    if (slots[tid] < 0)
      mooseError("Kokkos linear FV kernel '",
                 name(),
                 "': a diagonal matrix entry is missing from the sparsity pattern.");
  }
  _diagonal_matrix_slots = slots;
}

void
LinearFVElementalKernel::computeRightHandSide()
{
  Policy policy(0, numKokkosBlockElements());
  if (!_rhs_dispatcher)
    _rhs_dispatcher = DispatcherRegistry::build<RightHandSideLoop>(this, type());
  _rhs_dispatcher->parallelForAsync(policy);
}

void
LinearFVElementalKernel::computeMatrix()
{
  if (DispatcherRegistry::hasUserMethod<MatrixLoop>(type()))
  {
    Policy policy(0, numKokkosBlockElements());
    if (!_matrix_dispatcher)
      _matrix_dispatcher = DispatcherRegistry::build<MatrixLoop>(this, type());
    _matrix_dispatcher->parallelForAsync(policy);
  }
}

} // namespace Moose::Kokkos
