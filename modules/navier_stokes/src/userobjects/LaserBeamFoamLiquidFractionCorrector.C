//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LaserBeamFoamLiquidFractionCorrector.h"

#include "MooseMesh.h"
#include "MooseVariableFV.h"
#include "SystemBase.h"

#include <algorithm>
#include <cmath>

registerMooseObject("NavierStokesApp", LaserBeamFoamLiquidFractionCorrector);
registerMooseObjectAliased("NavierStokesApp",
                           LaserBeamFoamLiquidFractionCorrector,
                           "EnthalpyLiquidFractionCorrector");

InputParameters
LaserBeamFoamLiquidFractionCorrector::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params += BlockRestrictable::validParams();
  params.addClassDescription(
      "Applies an enthalpy liquid-fraction correction "
      "fl <- clamp(fl + relaxation * cp / L * (T - ((T_liquidus - T_solidus) * fl + "
      "T_solidus)), 0, 1).");
  params.addRequiredParam<VariableName>("liquid_fraction_variable",
                                        "Finite-volume liquid-fraction variable to update.");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature functor.");
  params.addRequiredParam<MooseFunctorName>("cp", "Specific heat functor.");
  params.addRequiredParam<MooseFunctorName>("L", "Latent heat functor.");
  params.addRequiredParam<MooseFunctorName>("T_solidus", "Solidus temperature.");
  params.addRequiredParam<MooseFunctorName>("T_liquidus", "Liquidus temperature.");
  params.addParam<MooseFunctorName>(
      "material_fraction",
      "Optional material volume fraction. When supplied, cells below min_material_fraction are not "
      "liquid-fraction corrected.");
  params.addRangeCheckedParam<Real>(
      "min_material_fraction",
      0.0,
      "min_material_fraction>=0 & min_material_fraction<=1",
      "Minimum material fraction required for liquid-fraction correction.");
  params.addRangeCheckedParam<Real>(
      "relaxation", 0.9, "relaxation>0 & relaxation<=1", "Liquid-fraction correction relaxation.");
  params.addRangeCheckedParam<Real>(
      "max_liquid_fraction_change",
      1.0,
      "max_liquid_fraction_change>0 & max_liquid_fraction_change<=1",
      "Maximum absolute time-step liquid-fraction change used to build the stored phase-change "
      "rate consumed by pressure sources. The thermodynamic liquid-fraction correction itself is "
      "not clipped by this parameter.");
  params.addRangeCheckedParam<Real>(
      "phase_change_rate_relaxation",
      1.0,
      "phase_change_rate_relaxation>0 & phase_change_rate_relaxation<=1",
      "Relaxation applied to the stored local phase-change rate consumed by pressure sources.");
  return params;
}

LaserBeamFoamLiquidFractionCorrector::LaserBeamFoamLiquidFractionCorrector(
    const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    BlockRestrictable(this),
    _liquid_fraction_var(dynamic_cast<MooseVariableFVReal &>(UserObject::_subproblem.getVariable(
        0, getParam<VariableName>("liquid_fraction_variable")))),
    _liquid_fraction_system(_liquid_fraction_var.sys()),
    _liquid_fraction_sys_num(_liquid_fraction_var.sys().number()),
    _liquid_fraction_var_num(_liquid_fraction_var.number()),
    _temperature(getFunctor<Real>("temperature")),
    _specific_heat(getFunctor<Real>("cp")),
    _latent_heat(getFunctor<Real>("L")),
    _T_solidus(getFunctor<Real>("T_solidus")),
    _T_liquidus(getFunctor<Real>("T_liquidus")),
    _material_fraction(isParamValid("material_fraction") ? &getFunctor<Real>("material_fraction")
                                                         : nullptr),
    _min_material_fraction(getParam<Real>("min_material_fraction")),
    _relaxation(getParam<Real>("relaxation")),
    _max_liquid_fraction_change(getParam<Real>("max_liquid_fraction_change")),
    _phase_change_rate_relaxation(getParam<Real>("phase_change_rate_relaxation"))
{
}

void
LaserBeamFoamLiquidFractionCorrector::execute()
{
  correctLiquidFraction();
}

Real
LaserBeamFoamLiquidFractionCorrector::correctLiquidFraction()
{
  auto & current_solution = *_liquid_fraction_system.system().current_local_solution;
  Real local_max = 0.0;
  Real local_weighted_sum = 0.0;
  Real local_volume = 0.0;

  const auto state = Moose::currentState();

  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!hasBlocks(elem_info->subdomain_id()))
      continue;

    const dof_id_type dof =
        elem_info->dofIndices()[_liquid_fraction_sys_num][_liquid_fraction_var_num];
    if (dof == DofObject::invalid_id)
      continue;

    const auto elem_arg = makeElemArg(elem_info->elem());
    if (_material_fraction && (*_material_fraction)(elem_arg, state) <= _min_material_fraction)
    {
      _phase_change_rate[dof] = 0.0;
      continue;
    }

    const Real latent_heat = _latent_heat(elem_arg, state);
    const Real delta_T = _T_liquidus(elem_arg, state) - _T_solidus(elem_arg, state);

    if (latent_heat <= 0.0 || delta_T <= 0.0)
    {
      _phase_change_rate[dof] = 0.0;
      continue;
    }

    const Real old_fl = std::clamp(current_solution(dof), 0.0, 1.0);
    const Real T_corr = delta_T * old_fl + _T_solidus(elem_arg, state);
    const Real proposed_fl = old_fl + _relaxation * _specific_heat(elem_arg, state) / latent_heat *
                                          (_temperature(elem_arg, state) - T_corr);
    const Real new_fl = std::clamp(proposed_fl, 0.0, 1.0);
    const Real correction = std::abs(new_fl - old_fl);
    const Real volume = elem_info->volume() * elem_info->coordFactor();

    current_solution.set(dof, new_fl);
    _liquid_fraction_system.solution().set(dof, new_fl);

    const Real old_time_fl = std::clamp(_liquid_fraction_system.solutionOld()(dof), 0.0, 1.0);
    const Real limited_rate_delta =
        std::clamp(new_fl - old_time_fl, -_max_liquid_fraction_change, _max_liquid_fraction_change);
    const Real unrelaxed_rate = _dt > 0.0 ? limited_rate_delta / _dt : 0.0;
    const auto rate_it = _phase_change_rate.find(dof);
    const Real old_rate = rate_it != _phase_change_rate.end() ? rate_it->second : 0.0;
    _phase_change_rate[dof] = (1.0 - _phase_change_rate_relaxation) * old_rate +
                              _phase_change_rate_relaxation * unrelaxed_rate;

    local_max = std::max(local_max, correction);
    local_weighted_sum += correction * volume;
    local_volume += volume;
  }

  current_solution.close();
  _liquid_fraction_system.solution().close();

  _fe_problem.comm().max(local_max);
  _fe_problem.comm().sum(local_weighted_sum);
  _fe_problem.comm().sum(local_volume);

  _max_correction = local_max;
  _mean_correction = local_volume > 0.0 ? local_weighted_sum / local_volume : 0.0;

  return _max_correction;
}

Real
LaserBeamFoamLiquidFractionCorrector::phaseChangeRate(const Elem * elem) const
{
  if (!elem)
    return 0.0;

  const auto & elem_info = _fe_problem.mesh().elemInfo(elem->id());
  const dof_id_type dof =
      elem_info.dofIndices()[_liquid_fraction_sys_num][_liquid_fraction_var_num];
  const auto rate_it = _phase_change_rate.find(dof);
  return rate_it != _phase_change_rate.end() ? rate_it->second : 0.0;
}

Real
LaserBeamFoamLiquidFractionCorrector::phaseChangeRate(const Moose::ElemArg & elem_arg) const
{
  return phaseChangeRate(elem_arg.elem);
}

Real
LaserBeamFoamLiquidFractionCorrector::phaseChangeRate(const Moose::FaceArg & face_arg) const
{
  if (face_arg.face_side)
    return phaseChangeRate(face_arg.face_side);

  const Real elem_rate = phaseChangeRate(face_arg.makeElem());
  if (!face_arg.fi || !face_arg.fi->neighborPtr())
    return elem_rate;

  return 0.5 * (elem_rate + phaseChangeRate(face_arg.makeNeighbor()));
}
