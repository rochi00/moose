//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LowMachEnthalpyNewtonState.h"

#include "LowMachEnthalpyThermodynamics.h"
#include "MooseLinearVariableFV.h"
#include "SystemBase.h"

#include <algorithm>
#include <cmath>

registerMooseObject("NavierStokesApp", LowMachEnthalpyNewtonState);

InputParameters
LowMachEnthalpyNewtonState::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params += BlockRestrictable::validParams();
  params.addClassDescription(
      "Restores the exact enthalpy-temperature relation after each Newton-linearized low-Mach "
      "temperature solve and measures convergence with the liquid-fraction L2 update.");
  params.addRequiredParam<VariableName>("temperature_variable",
                                        "Linear finite-volume temperature variable.");
  params.addRequiredParam<MooseFunctorName>("material_fraction", "Gas/PCM material fraction.");
  params.addRequiredParam<MooseFunctorName>("specific_enthalpy",
                                            "Exact specific enthalpy at the current temperature.");
  params.addRequiredParam<MooseFunctorName>("dh_dT", "Current Newton derivative of enthalpy.");
  params.addRequiredParam<MooseFunctorName>("liquid_fraction",
                                            "Exact current liquid-fraction functor.");
  params.addRequiredParam<MooseFunctorName>("cp_gas", "Gas specific heat.");
  params.addRequiredParam<MooseFunctorName>("cp_solid", "Solid PCM specific heat.");
  params.addRequiredParam<MooseFunctorName>("cp_liquid", "Liquid PCM specific heat.");
  params.addRequiredParam<MooseFunctorName>("rho_solid", "Solid PCM density.");
  params.addRequiredParam<MooseFunctorName>("rho_liquid", "Liquid PCM density.");
  params.addRequiredParam<MooseFunctorName>("h_solid", "Solidus specific enthalpy.");
  params.addRequiredParam<MooseFunctorName>("h_liquid", "Liquidus specific enthalpy.");
  params.addRequiredParam<MooseFunctorName>("T_solidus", "Solidus temperature.");
  params.addRequiredParam<MooseFunctorName>("T_liquidus", "Liquidus temperature.");
  params.addParam<MooseFunctorName>("T_reference", 0.0, "Reference temperature.");

  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_NONE);
  exec_enum = {EXEC_NONE};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

LowMachEnthalpyNewtonState::LowMachEnthalpyNewtonState(const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    BlockRestrictable(this),
    _temperature_variable(dynamic_cast<MooseLinearVariableFVReal &>(
        UserObject::_subproblem.getVariable(0, getParam<VariableName>("temperature_variable")))),
    _temperature_system(_temperature_variable.sys()),
    _system_name(_temperature_system.name()),
    _system_number(_temperature_system.number()),
    _variable_number(_temperature_variable.number()),
    _material_fraction_name(getParam<MooseFunctorName>("material_fraction")),
    _material_fraction(getFunctor<Real>("material_fraction")),
    _specific_enthalpy(getFunctor<Real>("specific_enthalpy")),
    _dh_dT(getFunctor<Real>("dh_dT")),
    _liquid_fraction(getFunctor<Real>("liquid_fraction")),
    _cp_gas(getFunctor<Real>("cp_gas")),
    _cp_solid(getFunctor<Real>("cp_solid")),
    _cp_liquid(getFunctor<Real>("cp_liquid")),
    _rho_solid(getFunctor<Real>("rho_solid")),
    _rho_liquid(getFunctor<Real>("rho_liquid")),
    _h_solid(getFunctor<Real>("h_solid")),
    _h_liquid(getFunctor<Real>("h_liquid")),
    _T_solidus(getFunctor<Real>("T_solidus")),
    _T_liquidus(getFunctor<Real>("T_liquidus")),
    _T_reference(getFunctor<Real>("T_reference"))
{
}

void
LowMachEnthalpyNewtonState::beginNewtonIteration()
{
  _iteration_state.clear();
  _has_captured_update = false;
  const auto state = Moose::currentState();
  const auto & solution = _temperature_system.solution();

  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!hasBlocks(elem_info->subdomain_id()))
      continue;

    const dof_id_type dof = elem_info->dofIndices()[_system_number][_variable_number];
    if (dof == DofObject::invalid_id || dof < solution.first_local_index() ||
        dof >= solution.last_local_index())
      continue;

    const auto elem_arg = makeElemArg(elem_info->elem());
    _iteration_state.emplace(dof,
                             IterationState{solution(dof),
                                            _specific_enthalpy(elem_arg, state),
                                            _dh_dT(elem_arg, state),
                                            _liquid_fraction(elem_arg, state),
                                            _material_fraction(elem_arg, state),
                                            _cp_gas(elem_arg, state),
                                            _cp_solid(elem_arg, state),
                                            _cp_liquid(elem_arg, state),
                                            _rho_solid(elem_arg, state),
                                            _rho_liquid(elem_arg, state),
                                            _h_solid(elem_arg, state),
                                            _h_liquid(elem_arg, state),
                                            _T_solidus(elem_arg, state),
                                            _T_liquidus(elem_arg, state),
                                            _T_reference(elem_arg, state),
                                            0.0});
  }
}

void
LowMachEnthalpyNewtonState::captureNewtonUpdate()
{
  const auto & solution = _temperature_system.solution();

  for (auto & [dof, old_state] : _iteration_state)
  {
    const Real linear_temperature = solution(dof);
    old_state.full_step_enthalpy = LowMachEnthalpy::newtonUpdatedEnthalpy(
        old_state.enthalpy, old_state.dh_dT, old_state.temperature, linear_temperature);

    if (!std::isfinite(old_state.full_step_enthalpy))
      mooseError(name(), ": non-finite Newton enthalpy produced at temperature DOF ", dof, ".");
  }

  _has_captured_update = true;
}

Real
LowMachEnthalpyNewtonState::restoreExactState(const Real step_length)
{
  if (!_has_captured_update)
    captureNewtonUpdate();
  if (step_length <= 0.0 || step_length > 1.0)
    mooseError(name(), ": Newton step length must be in (0, 1].");

  auto & solution = _temperature_system.solution();
  auto & current_local_solution = *_temperature_system.system().current_local_solution;
  Real local_delta_norm_squared = 0.0;
  Real local_old_norm_squared = 0.0;

  for (const auto & [dof, old_state] : _iteration_state)
  {
    const Real updated_enthalpy =
        old_state.enthalpy + step_length * (old_state.full_step_enthalpy - old_state.enthalpy);
    const Real updated_temperature = exactTemperature(updated_enthalpy, old_state);
    const Real updated_liquid_fraction = exactLiquidFraction(updated_enthalpy, old_state);

    if (!std::isfinite(updated_enthalpy) || !std::isfinite(updated_temperature) ||
        !std::isfinite(updated_liquid_fraction))
      mooseError(name(), ": non-finite thermodynamic state produced at temperature DOF ", dof, ".");

    solution.set(dof, updated_temperature);
    current_local_solution.set(dof, updated_temperature);

    const Real liquid_fraction_update = updated_liquid_fraction - old_state.liquid_fraction;
    local_delta_norm_squared += liquid_fraction_update * liquid_fraction_update;
    local_old_norm_squared += old_state.liquid_fraction * old_state.liquid_fraction;
  }

  solution.close();
  current_local_solution.close();

  _communicator.sum(local_delta_norm_squared);
  _communicator.sum(local_old_norm_squared);
  return std::sqrt(local_delta_norm_squared) / (1.0 + std::sqrt(local_old_norm_squared));
}

Real
LowMachEnthalpyNewtonState::exactTemperature(const Real enthalpy,
                                             const IterationState & state) const
{
  if (state.material_fraction < 0.5)
    return enthalpy / state.cp_gas + state.T_reference;
  if (enthalpy < state.h_solid)
    return enthalpy / state.cp_solid + state.T_reference;
  if (enthalpy <= state.h_liquid)
    return state.T_solidus + (enthalpy - state.h_solid) * (state.T_liquidus - state.T_solidus) /
                                 (state.h_liquid - state.h_solid);
  return state.T_liquidus + (enthalpy - state.h_liquid) / state.cp_liquid;
}

Real
LowMachEnthalpyNewtonState::exactLiquidFraction(const Real enthalpy,
                                                const IterationState & state) const
{
  if (state.material_fraction < 0.5 || enthalpy < state.h_solid)
    return 0.0;
  if (enthalpy > state.h_liquid)
    return 1.0;

  const Real denominator = enthalpy * (state.rho_liquid - state.rho_solid) -
                           state.rho_liquid * state.h_liquid + state.rho_solid * state.h_solid;
  return std::clamp(state.rho_solid * (state.h_solid - enthalpy) / denominator, 0.0, 1.0);
}
