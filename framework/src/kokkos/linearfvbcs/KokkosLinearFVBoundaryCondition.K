//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVBoundaryCondition.h"

#include "FEProblemBase.h"
#include "InputParameters.h"
#include "SubProblem.h"
#include "SystemBase.h"

namespace Moose::Kokkos
{

InputParameters
LinearFVBoundaryCondition::validParams()
{
  auto params = MooseObject::validParams();
  params += TransientInterface::validParams();
  params += MeshChangedInterface::validParams();
  params += BoundaryRestrictableRequired::validParams();

  params.addPrivateParam<bool>(MooseBase::kokkos_object_param, true);
  // The objects read the cell values of variables on the one-layer ghost cells, including those of
  // auxiliary variables, whose system has no matrix coupling that would ghost them
  params.addRelationshipManager("ElementSideNeighborLayers",
                                Moose::RelationshipManagerType::ALGEBRAIC,
                                [](const InputParameters &, InputParameters & rm_params)
                                { rm_params.set<unsigned short>("layers") = 1; });
  params.addRequiredParam<LinearVariableName>(
      "variable", "The name of the variable this boundary condition supplies data for");
  params.declareControllable("enable");
  params.set<bool>("_residual_object") = false;
  params.registerBase("LinearFVBoundaryCondition");
  params.registerSystemAttributeName("LinearFVBoundaryCondition");
  return params;
}

LinearFVBoundaryCondition::LinearFVBoundaryCondition(const InputParameters & parameters)
  : MooseObject(parameters),
    SetupInterface(this),
    FunctionInterface(this),
    UserObjectInterface(this),
    TransientInterface(this),
    PostprocessorInterface(this),
    VectorPostprocessorInterface(this, false),
    Restartable(this, getBase() + "s"),
    MeshChangedInterface(parameters),
    BoundaryRestrictableRequired(this, false),
    MeshHolder(*getCheckedPointerParam<SubProblem *>("_subproblem")->mesh().getKokkosMesh()),
    FVMeshHolder(getCheckedPointerParam<SystemBase *>("_sys")->feProblem().getKokkosFVMesh()),
    SystemHolder(getCheckedPointerParam<SystemBase *>("_sys")->feProblem().getKokkosSystems()),
    LinearFVVariableInterface(*this,
                              *getParam<FEProblemBase *>("_fe_problem_base"),
                              getCheckedPointerParam<SystemBase *>("_sys")->number()),
    LinearFVPropertyInterface(*this, *getParam<FEProblemBase *>("_fe_problem_base"), *this),
    _fe_problem(*getParam<FEProblemBase *>("_fe_problem_base")),
    _t(TransientInterface::_t),
    _t_old(TransientInterface::_t_old),
    _t_step(TransientInterface::_t_step),
    _dt(TransientInterface::_dt),
    _dt_old(TransientInterface::_dt_old)
{
  const auto * const subproblem = getCheckedPointerParam<SubProblem *>("_subproblem");
  _var.init(subproblem->getVariable(0,
                                    getParam<LinearVariableName>("variable"),
                                    Moose::VarKindType::VAR_SOLVER,
                                    Moose::VarFieldType::VAR_FIELD_STANDARD),
            Moose::SOLUTION_TAG);
}

LinearFVBoundaryCondition::LinearFVBoundaryCondition(const LinearFVBoundaryCondition & object)
  : MooseObject(object, {}),
    SetupInterface(object, {}),
    FunctionInterface(object, {}),
    UserObjectInterface(object, {}),
    TransientInterface(object, {}),
    PostprocessorInterface(object, {}),
    VectorPostprocessorInterface(object, {}),
    Restartable(object, {}),
    MeshChangedInterface(object, {}),
    BoundaryRestrictableRequired(object, {}),
    MeshHolder(object),
    FVMeshHolder(object),
    SystemHolder(object),
    LinearFVVariableInterface(object),
    LinearFVPropertyInterface(object),
    _fe_problem(object._fe_problem),
    _var(object._var),
    _boundary_data(object._boundary_data),
    _faces(object._faces),
    _face_range_flips(object._face_range_flips),
    _t(object._t),
    _t_old(object._t_old),
    _t_step(object._t_step),
    _dt(object._dt),
    _dt_old(object._dt_old)
{
}

void
LinearFVBoundaryCondition::initialSetup()
{
  setupFaces();
  setupLinearFVProperties();
}

void
LinearFVBoundaryCondition::setupFaces()
{
  if (_faces_built)
    return;
  _faces_built = true;

  const auto & fv_mesh = kokkosFVMesh();
  const auto * const variable = _var.mooseFieldVar();
  const auto & bc_boundaries = boundaryIDs();

  std::vector<FaceRange> ranges;
  std::vector<bool> flips;
  for (const auto r : index_range(fv_mesh.getFaceRanges()))
  {
    const auto & range = fv_mesh.getFaceRanges()[r];
    const bool on_elem = variable->activeOnSubdomain(range.elem_subdomain);
    const bool on_neighbor = range.neighbor_subdomain != Moose::INVALID_BLOCK_ID &&
                             variable->activeOnSubdomain(range.neighbor_subdomain);
    // A boundary for the variable, on the side of an owned cell
    if (on_elem == on_neighbor || !(on_elem ? range.elem_owned : range.neighbor_owned))
      continue;

    // The boundary IDs that select the variable's boundary condition
    const auto & range_boundaries = fv_mesh.getFaceRangeVariableBoundaryIDs(r, on_neighbor);
    if (range_boundaries.size() != 1 || !bc_boundaries.count(range_boundaries.front()))
      continue;

    _mesh_range_to_range[r] = ranges.size();
    ranges.push_back(range);
    flips.push_back(on_neighbor);
  }

  _faces.create(ranges);
  _host_face_range_flips = flips;
  _face_range_flips.create(flips.size());
  for (const auto i : index_range(flips))
    _face_range_flips[i] = flips[i];
  _face_range_flips.copyToDevice();

  const auto n_boundary_faces = numBoundaryFaces();

  _boundary_data.value.coefficient.create(n_boundary_faces);
  _boundary_data.value.source.create(n_boundary_faces);
  _boundary_data.normal_gradient.coefficient.create(n_boundary_faces);
  _boundary_data.normal_gradient.source.create(n_boundary_faces);
  _boundary_data.value.coefficient = 0;
  _boundary_data.value.source = 0;
  _boundary_data.normal_gradient.coefficient = 0;
  _boundary_data.normal_gradient.source = 0;
}

std::optional<MOOSE_KOKKOS_INDEX_TYPE>
LinearFVBoundaryCondition::faceRangeOffset(const unsigned int range) const
{
  mooseAssert(_faces_built, "setupFaces() must be called first");

  const auto it = _mesh_range_to_range.find(range);
  if (it == _mesh_range_to_range.end())
    return {};
  return _faces.offset(it->second);
}

bool
LinearFVBoundaryCondition::hasBoundaryValue() const
{
  return DispatcherRegistry::hasUserMethod<BoundaryValueLoop>(type());
}

bool
LinearFVBoundaryCondition::hasBoundaryNormalGradient() const
{
  return DispatcherRegistry::hasUserMethod<BoundaryNormalGradientLoop>(type());
}

void
LinearFVBoundaryCondition::computeBoundaryValueData()
{
  if (!hasBoundaryValue())
    mooseError("Kokkos linear FV boundary condition '",
               name(),
               "' (type '",
               type(),
               "') was asked to compute boundary values, but it does not override "
               "computeBoundaryValue().");

  updateLinearFVProperties(_t, nullptr, &_faces, _host_face_range_flips);
  Policy policy(0, numBoundaryFaces());
  if (!_boundary_value_dispatcher)
    _boundary_value_dispatcher = DispatcherRegistry::build<BoundaryValueLoop>(this, type());
  _boundary_value_dispatcher->parallelForAsync(policy);
}

void
LinearFVBoundaryCondition::computeBoundaryNormalGradientData()
{
  if (!hasBoundaryNormalGradient())
    mooseError("Kokkos linear FV boundary condition '",
               name(),
               "' (type '",
               type(),
               "') was asked to compute boundary normal gradients, but it does not override "
               "computeBoundaryNormalGradient().");

  updateLinearFVProperties(_t, nullptr, &_faces, _host_face_range_flips);
  Policy policy(0, numBoundaryFaces());
  if (!_boundary_normal_gradient_dispatcher)
    _boundary_normal_gradient_dispatcher =
        DispatcherRegistry::build<BoundaryNormalGradientLoop>(this, type());
  _boundary_normal_gradient_dispatcher->parallelForAsync(policy);
}

} // namespace Moose::Kokkos
