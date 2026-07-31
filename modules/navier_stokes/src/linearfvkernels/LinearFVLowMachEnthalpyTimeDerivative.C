//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVLowMachEnthalpyTimeDerivative.h"

#include "MooseLinearVariableFV.h"

registerMooseObject("NavierStokesApp", LinearFVLowMachEnthalpyTimeDerivative);

InputParameters
LinearFVLowMachEnthalpyTimeDerivative::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Adds the Newton-linearized conservative low-Mach enthalpy time derivative to a "
      "temperature equation.");
  params.addRequiredParam<MooseFunctorName>(
      "temporary_density", "Fixed end-of-step density produced by the conservative mass solve.");
  params.addRequiredParam<MooseFunctorName>(
      "density", "Accepted EOS density; its old state supplies the old conservative density.");
  params.addRequiredParam<MooseFunctorName>("specific_enthalpy",
                                            "Exact specific enthalpy evaluated from temperature.");
  params.addRequiredParam<MooseFunctorName>("dh_dT", "Current Newton derivative of enthalpy.");
  return params;
}

LinearFVLowMachEnthalpyTimeDerivative::LinearFVLowMachEnthalpyTimeDerivative(
    const InputParameters & params)
  : LinearFVElementalKernel(params),
    _temporary_density(getFunctor<Real>("temporary_density")),
    _density(getFunctor<Real>("density")),
    _specific_enthalpy(getFunctor<Real>("specific_enthalpy")),
    _dh_dT(getFunctor<Real>("dh_dT")),
    _time_integrator(_sys.getTimeIntegrator(_var_num))
{
}

Real
LinearFVLowMachEnthalpyTimeDerivative::computeMatrixContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();
  const auto coefficients = _time_integrator.timeDerivativeCoefficients();
  return coefficients.front() * _temporary_density(elem_arg, state) * _dh_dT(elem_arg, state) /
         _dt * _current_elem_volume;
}

Real
LinearFVLowMachEnthalpyTimeDerivative::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();

  const Real temporary_density = _temporary_density(elem_arg, state);
  const Real enthalpy = _specific_enthalpy(elem_arg, state);
  const Real dh_dT = _dh_dT(elem_arg, state);
  const Real temperature = _var.getElemValue(*_current_elem_info, state);
  const auto coefficients = _time_integrator.timeDerivativeCoefficients();

  Real temporal_residual =
      coefficients.front() * temporary_density * (enthalpy - dh_dT * temperature);
  for (const auto state_index : index_range(coefficients))
  {
    if (state_index == 0)
      continue;

    const Moose::StateArg history_state(state_index, Moose::SolutionIterationType::Time);
    temporal_residual += coefficients[state_index] * _density(elem_arg, history_state) *
                         _specific_enthalpy(elem_arg, history_state);
  }

  return -temporal_residual / _dt * _current_elem_volume;
}
