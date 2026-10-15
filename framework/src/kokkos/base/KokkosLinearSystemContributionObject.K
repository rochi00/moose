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

#include <algorithm>

namespace Moose::Kokkos
{

InputParameters
LinearSystemContributionObject::validParams()
{
  auto params = ::LinearSystemContributionObject::validParams();
  params.addPrivateParam<bool>(MooseBase::kokkos_object_param, true);
  // The objects read the cell values of variables on the one-layer ghost cells, including those of
  // auxiliary variables, whose system has no matrix coupling that would ghost them
  params.addRelationshipManager("ElementSideNeighborLayers",
                                Moose::RelationshipManagerType::ALGEBRAIC,
                                [](const InputParameters &, InputParameters & rm_params)
                                { rm_params.set<unsigned short>("layers") = 1; });
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
    FVMeshHolder(getCheckedPointerParam<SystemBase *>("_sys")->feProblem().getKokkosFVMesh()),
    SystemHolder(getCheckedPointerParam<SystemBase *>("_sys")->feProblem().getKokkosSystems()),
    LinearFVVariableInterface(*this,
                              *getParam<FEProblemBase *>("_fe_problem_base"),
                              getCheckedPointerParam<SystemBase *>("_sys")->number()),
    LinearFVPropertyInterface(*this, *getParam<FEProblemBase *>("_fe_problem_base"), *this),
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
    FVMeshHolder(object),
    SystemHolder(object),
    LinearFVVariableInterface(object),
    LinearFVPropertyInterface(object),
    _fe_problem(object._fe_problem),
    _var(object._var),
    _assemble_matrix(object._assemble_matrix),
    _assemble_rhs(object._assemble_rhs),
    _num_vector_targets(object._num_vector_targets),
    _num_matrix_targets(object._num_matrix_targets),
    _t(object._t)
{
  std::copy_n(object._vector_targets, max_tagged_targets, _vector_targets);
  std::copy_n(object._matrix_targets, max_tagged_targets, _matrix_targets);
}

void
LinearSystemContributionObject::computeContributions(const bool matrix,
                                                     const bool rhs,
                                                     const FVSystemAssembly & assembly)
{
  // The arrays of this object's tags that are being assembled. Tags can be added to the governing
  // TaggingInterface sets after construction, so they are read at every assembly.
  auto collect = [this](const std::set<TagID> & tags,
                        const bool assemble,
                        auto && values,
                        Real ** targets,
                        unsigned int & num_targets)
  {
    num_targets = 0;
    if (!assemble)
      return;
    for (const auto tag : tags)
      if (auto * const array = values(tag))
      {
        if (num_targets == max_tagged_targets)
          mooseError("Kokkos linear FV object '",
                     name(),
                     "' adds to more than ",
                     max_tagged_targets,
                     " tagged matrices or vectors in one assembly, which is not supported.");
        targets[num_targets++] = array;
      }
  };
  collect(getMatrixTags({}),
          matrix,
          [&](const TagID tag) { return assembly.matrixValues(tag); },
          _matrix_targets,
          _num_matrix_targets);
  collect(getVectorTags({}),
          rhs,
          [&](const TagID tag) { return assembly.vectorValues(tag); },
          _vector_targets,
          _num_vector_targets);

  _assemble_matrix = _num_matrix_targets;
  _assemble_rhs = _num_vector_targets;
  if (_assemble_matrix || _assemble_rhs)
    dispatchContributions();
}

} // namespace Moose::Kokkos
