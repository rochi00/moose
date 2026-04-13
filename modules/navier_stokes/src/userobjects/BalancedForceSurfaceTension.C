//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BalancedForceSurfaceTension.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "ElemInfo.h"
#include "FaceInfo.h"

#include "libmesh/linear_implicit_system.h"

registerMooseObject("NavierStokesApp", BalancedForceSurfaceTension);

InputParameters
BalancedForceSurfaceTension::validParams()
{
  auto params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params.addClassDescription(
      "Provides curvature and alpha values for balanced-force surface tension "
      "in the Rhie-Chow face flux. The curvature is read from an existing "
      "AuxVariable (e.g., computed by LinearFVCurvatureAux) and interpolated "
      "to faces. Garbage curvature far from the interface is harmless because "
      "it is multiplied by (alpha_N - alpha_C) = 0 in the face flux.");
  params.addRequiredParam<Real>("sigma", "Surface tension coefficient.");
  params.addRequiredParam<MooseFunctorName>("kappa", "Curvature field (cell-centered functor).");
  params.addRequiredParam<VariableName>("alpha_variable", "The VOF phase fraction variable.");
  return params;
}

BalancedForceSurfaceTension::BalancedForceSurfaceTension(const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    _sigma(getParam<Real>("sigma")),
    _kappa(getFunctor<Real>("kappa")),
    _alpha_var(dynamic_cast<MooseLinearVariableFV<Real> &>(
        _fe_problem.getVariable(0, getParam<VariableName>("alpha_variable"))))
{
}

Real
BalancedForceSurfaceTension::getFaceCurvature(const FaceInfo & fi) const
{
  const Moose::FaceArg face_arg{
      &fi, Moose::FV::LimiterType::CentralDifference, true, false, fi.elemPtr(), nullptr};
  return MetaPhysicL::raw_value(_kappa(face_arg, Moose::currentState()));
}

Real
BalancedForceSurfaceTension::getAlpha(dof_id_type elem_id) const
{
  const auto & alpha_sys =
      libMesh::cast_ref<const libMesh::LinearImplicitSystem &>(_alpha_var.sys().system());
  const auto & ei = _fe_problem.mesh().elemInfo(elem_id);
  const auto alpha_dof = ei.dofIndices()[_alpha_var.sys().number()][_alpha_var.number()];
  return (*alpha_sys.current_local_solution)(alpha_dof);
}
