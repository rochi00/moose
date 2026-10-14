//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearSystemContributionObject.h"

#include "FEProblemBase.h"
#include "InputParameters.h"

namespace Moose::Kokkos
{

InputParameters
LinearSystemContributionObject::validParams()
{
  auto params = ::LinearSystemContributionObject::validParams();
  params.addPrivateParam<bool>(MooseBase::kokkos_object_param, true);
  return params;
}

LinearSystemContributionObject::LinearSystemContributionObject(const InputParameters & parameters)
  : MooseObject(parameters),
    SetupInterface(this),
    FunctionInterface(this),
    UserObjectInterface(this),
    TransientInterface(this),
    PostprocessorInterface(this),
    VectorPostprocessorInterface(this, false),
    Restartable(this, getBase() + "s"),
    MeshChangedInterface(parameters),
    TaggingInterface(this),
    MeshHolder(*getCheckedPointerParam<SubProblem *>("_subproblem")->mesh().getKokkosMesh()),
    SystemHolder(getCheckedPointerParam<SystemBase *>("_sys")->feProblem().getKokkosSystems()),
    _fe_problem(*getParam<FEProblemBase *>("_fe_problem_base")),
    _t(TransientInterface::_t)
{
  const auto * const subproblem = getCheckedPointerParam<SubProblem *>("_subproblem");
  _var.init(subproblem->getVariable(0,
                                    getParam<LinearVariableName>("variable"),
                                    Moose::VarKindType::VAR_SOLVER,
                                    Moose::VarFieldType::VAR_FIELD_STANDARD),
            Moose::SOLUTION_TAG);
}

LinearSystemContributionObject::LinearSystemContributionObject(
    const LinearSystemContributionObject & object)
  : MooseObject(object, {}),
    SetupInterface(object, {}),
    FunctionInterface(object, {}),
    UserObjectInterface(object, {}),
    TransientInterface(object, {}),
    PostprocessorInterface(object, {}),
    VectorPostprocessorInterface(object, {}),
    Restartable(object, {}),
    MeshChangedInterface(object, {}),
    TaggingInterface(object, {}),
    MeshHolder(object),
    SystemHolder(object),
    _fe_problem(object._fe_problem),
    _var(object._var),
    _t(object._t)
{
  // Tags may be added to the governing tag sets held by the TaggingInterface later than the time
  // that the primary constructor is run. So it's more robust to build our tag set, immediately
  // before being dispatched to device, from the governing TaggingInterface set instead of relying
  // on the tag sets at construction
  object.updateTagArrays();
  _vector_tags = object._vector_tags;
  _matrix_tags = object._matrix_tags;
}

void
LinearSystemContributionObject::updateTagArrays() const
{
  const auto & vector_tags = getVectorTags({});
  const auto & matrix_tags = getMatrixTags({});

  if (_tag_arrays_built && vector_tags == _vector_tag_set && matrix_tags == _matrix_tag_set)
    return;

  _vector_tags = vector_tags;
  _matrix_tags = matrix_tags;
  _vector_tag_set = vector_tags;
  _matrix_tag_set = matrix_tags;
  _tag_arrays_built = true;
}

} // namespace Moose::Kokkos
