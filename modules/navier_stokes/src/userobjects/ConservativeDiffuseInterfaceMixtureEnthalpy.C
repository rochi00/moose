//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConservativeDiffuseInterfaceMixtureEnthalpy.h"

#include "SinglePhaseFluidProperties.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "ElemInfo.h"
#include "SystemBase.h"

#include <limits>

registerMooseObject("NavierStokesApp", ConservativeDiffuseInterfaceMixtureEnthalpy);

namespace
{
/// Relative temperature change at which the Newton recovery of T stops, near round-off for the
/// temperatures of order 1e2-1e3 K this is used at
constexpr Real temperature_relative_tolerance = 1e-12;
/// Newton iterations allowed per cell; h_k(T) is smooth and monotonic, so the recovery converges
/// quadratically in a few iterations from the lagged temperature. Where bisection takes over on a
/// melting interval, it reaches the tolerance from a bracket of order 1 K in about 32 halvings.
constexpr unsigned int max_newton_iterations = 50;
/// Temperature width over which the linearization alone spreads the latent heat of an isothermal
/// melt, so that dT/d(rho*h) of a cell at the melting temperature stays finite. Held exactly
/// there, a cell whose latent content is small against the heat it exchanges within a step, such
/// as a gas cell with a trace of liquid, gives up far more enthalpy than it holds and leaves the
/// recovery hundreds of kelvin below the solidus. At 0.01 K the effective heat capacity of a
/// liquid cell is still thousands of times its sensible one, and the recovery stays exact.
constexpr Real isothermal_linearization_width = 1e-2;
}

InputParameters
ConservativeDiffuseInterfaceMixtureEnthalpy::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params += BlockRestrictable::validParams();
  params.addClassDescription(
      "Relates the conserved mixture enthalpy density of a two-phase diffuse-interface flow to its "
      "temperature, sets it initially and recovers the temperature from it.");
  params.addRequiredParam<VariableName>("enthalpy_variable",
                                        "The solved mixture enthalpy density rho*h.");
  params.addRequiredParam<VariableName>("temperature",
                                        "The temperature auxiliary variable recovered from rho*h.");
  params.addRequiredParam<MooseFunctorName>("volume_fraction", "The liquid volume fraction.");
  params.addRequiredParam<MooseFunctorName>("liquid_density", "The constant liquid density.");
  params.addRequiredParam<MooseFunctorName>("gas_density", "The constant gas density.");
  params.addRequiredParam<UserObjectName>("liquid_fluid_properties",
                                          "Fluid properties giving the liquid enthalpy h(p, T).");
  params.addRequiredParam<UserObjectName>("gas_fluid_properties",
                                          "Fluid properties giving the gas enthalpy h(p, T).");
  params.addRequiredParam<Real>("thermodynamic_pressure",
                                "Pressure at which the phase enthalpies are evaluated.");

  params.addRangeCheckedParam<Real>("latent_heat",
                                    0,
                                    "latent_heat >= 0",
                                    "Latent heat of fusion of the liquid phase [J/kg]. The liquid "
                                    "phase melts and solidifies when it is positive.");
  params.addParam<Real>("solidus_temperature",
                        "Temperature below which the liquid phase is fully solid [K].");
  params.addParam<Real>("liquidus_temperature",
                        "Temperature above which the liquid phase is fully liquid [K]. Equal to "
                        "the solidus temperature for melting at a single temperature.");
  MooseEnum liquid_fraction_shape("linear smoothstep", "linear");
  params.addParam<MooseEnum>("liquid_fraction_shape",
                             liquid_fraction_shape,
                             "Liquid fraction over a finite melting interval: linear, or the "
                             "smoothstep 3 s^2 - 2 s^3 of the scaled temperature s.");
  params.addParam<VariableName>(
      "liquid_fraction",
      "The liquid fraction auxiliary variable recovered from rho*h, required with latent heat.");
  params.addParam<UserObjectName>(
      "solid_fluid_properties",
      "Fluid properties giving the heat capacity and thermal conductivity of the liquid phase's "
      "solid. The solid shares the liquid's own when omitted.");
  params.addParam<MooseFunctorName>(
      "solid_density",
      "The constant density of the liquid phase's solid. The solid shares the liquid density "
      "when omitted.");
  MooseEnum conductivity_mixing("series parallel", "series");
  params.addParam<MooseEnum>("conductivity_mixing",
                             conductivity_mixing,
                             "How the thermal conductivities of the phases combine in a cell: in "
                             "series, 1/k = sum_k alpha_k / k_k, or in parallel, k = sum_k alpha_k "
                             "k_k.");

  // Only the initial rho*h is set on execution; the solve recovers the temperature explicitly
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum = {EXEC_INITIAL};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  // Before the initial postprocessors, which may read rho*h
  params.set<int>("execution_order_group") = -1;
  return params;
}

ConservativeDiffuseInterfaceMixtureEnthalpy::ConservativeDiffuseInterfaceMixtureEnthalpy(
    const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    BlockRestrictable(this),
    _enthalpy_name(getParam<VariableName>("enthalpy_variable")),
    _temperature_name(getParam<VariableName>("temperature")),
    _alpha(getFunctor<Real>("volume_fraction")),
    _density({&getFunctor<Real>("liquid_density"), &getFunctor<Real>("gas_density")}),
    _solid_density(isParamValid("solid_density") ? &getFunctor<Real>("solid_density")
                                                 : _density[LIQUID]),
    _fp({&getUserObject<SinglePhaseFluidProperties>("liquid_fluid_properties"),
         &getUserObject<SinglePhaseFluidProperties>("gas_fluid_properties")}),
    _solid_fp(isParamValid("solid_fluid_properties")
                  ? &getUserObject<SinglePhaseFluidProperties>("solid_fluid_properties")
                  : _fp[LIQUID]),
    _thermodynamic_pressure(getParam<Real>("thermodynamic_pressure")),
    _latent_heat(getParam<Real>("latent_heat")),
    _solidus_temperature(isParamValid("solidus_temperature") ? getParam<Real>("solidus_temperature")
                                                             : 0.0),
    _liquidus_temperature(
        isParamValid("liquidus_temperature") ? getParam<Real>("liquidus_temperature") : 0.0),
    _smoothstep(getParam<MooseEnum>("liquid_fraction_shape") == "smoothstep"),
    _parallel_conductivity(getParam<MooseEnum>("conductivity_mixing") == "parallel")
{
  if (_latent_heat > 0)
  {
    for (const auto & param : {"solidus_temperature", "liquidus_temperature", "liquid_fraction"})
      if (!isParamValid(param))
        paramError(param, "Required when there is latent heat.");
    if (_liquidus_temperature < _solidus_temperature)
      paramError("liquidus_temperature",
                 "The liquidus temperature must not be below the solidus temperature.");
    // The two fluid properties objects have reference states of their own; the latent heat must be
    // the whole enthalpy difference between the phases at the solidus
    Real h_s, dh_dp, dh_dT;
    _solid_fp->h_from_p_T(_thermodynamic_pressure, _solidus_temperature, h_s, dh_dp, dh_dT);
    _solid_enthalpy_offset = phaseEnthalpy(LIQUID, _solidus_temperature).first - h_s;
  }
  else
    for (const auto & param : {"solidus_temperature",
                               "liquidus_temperature",
                               "liquid_fraction_shape",
                               "liquid_fraction",
                               "solid_fluid_properties",
                               "solid_density"})
      if (isParamSetByUser(param))
        paramError(param, "Only used when there is latent heat.");
}

void
ConservativeDiffuseInterfaceMixtureEnthalpy::initialSetup()
{
  cacheVariables();
}

void
ConservativeDiffuseInterfaceMixtureEnthalpy::meshChanged()
{
  cacheVariables();
}

void
ConservativeDiffuseInterfaceMixtureEnthalpy::cacheVariables()
{
  _enthalpy_var = dynamic_cast<MooseLinearVariableFVReal *>(
      &UserObject::_subproblem.getVariable(_tid, _enthalpy_name));
  if (!_enthalpy_var)
    paramError("enthalpy_variable", "The enthalpy variable must be a MooseLinearVariableFVReal.");

  _temperature_var = dynamic_cast<MooseLinearVariableFVReal *>(
      &UserObject::_subproblem.getVariable(_tid, _temperature_name));
  if (!_temperature_var)
    paramError("temperature", "The temperature variable must be a MooseLinearVariableFVReal.");

  if (isParamValid("liquid_fraction"))
  {
    _liquid_fraction_var = dynamic_cast<MooseLinearVariableFVReal *>(
        &UserObject::_subproblem.getVariable(_tid, getParam<VariableName>("liquid_fraction")));
    if (!_liquid_fraction_var)
      paramError("liquid_fraction",
                 "The liquid fraction variable must be a MooseLinearVariableFVReal.");
  }
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::cellAlpha(const ElemInfo & elem_info) const
{
  return _alpha(makeElemArg(elem_info.elem()), Moose::currentState());
}

ConservativeDiffuseInterfaceMixtureEnthalpy::Densities
ConservativeDiffuseInterfaceMixtureEnthalpy::cellDensities(const Moose::ElemArg & elem_arg,
                                                           const Moose::StateArg & state) const
{
  return {(*_density[LIQUID])(elem_arg, state),
          (*_solid_density)(elem_arg, state),
          (*_density[GAS])(elem_arg, state)};
}

std::pair<Real, Real>
ConservativeDiffuseInterfaceMixtureEnthalpy::phaseEnthalpy(const Phase phase,
                                                           const Real temperature) const
{
  Real h, dh_dp, dh_dT;
  _fp[phase]->h_from_p_T(_thermodynamic_pressure, temperature, h, dh_dp, dh_dT);
  return {h, dh_dT};
}

std::pair<Real, Real>
ConservativeDiffuseInterfaceMixtureEnthalpy::liquidFraction(const Real temperature) const
{
  mooseAssert(!isothermal(), "An isothermal melt has no liquid fraction function of temperature");
  const Real dT_melt = _liquidus_temperature - _solidus_temperature;
  const Real s = std::clamp((temperature - _solidus_temperature) / dT_melt, 0.0, 1.0);
  if (_smoothstep)
    return {s * s * (3.0 - 2.0 * s), 6.0 * s * (1.0 - s) / dT_melt};
  // The derivative of the linear fraction is a top hat over the open melting interval
  return {s, (s > 0.0 && s < 1.0) ? 1.0 / dT_melt : 0.0};
}

bool
ConservativeDiffuseInterfaceMixtureEnthalpy::hasPhaseChange(const Real liquid_partial_density,
                                                            const Real enthalpy_scale) const
{
  // rho*h carries the absolute enthalpy of the phases, so its round-off is of order
  // epsilon |rho*h|. The square root keeps the liquid fraction, the latent part of rho*h divided by
  // alpha rho_l L, accurate to about 1e-8 in the cells that change phase. Below that a cell holds
  // a negligible latent heat, and its liquid fraction is carried along unchanged.
  static const Real relative_resolution = std::sqrt(std::numeric_limits<Real>::epsilon());
  return liquid_partial_density * _latent_heat > relative_resolution * std::abs(enthalpy_scale);
}

ConservativeDiffuseInterfaceMixtureEnthalpy::EnthalpyWithDerivatives
ConservativeDiffuseInterfaceMixtureEnthalpy::liquidPhaseEnthalpyDensity(
    const Real temperature, const Real liquid_fraction, const Densities & densities) const
{
  const auto [h_l, cp_l] = phaseEnthalpy(LIQUID, temperature);
  Real h_s, dh_dp, cp_s;
  _solid_fp->h_from_p_T(_thermodynamic_pressure, temperature, h_s, dh_dp, cp_s);
  h_s += _solid_enthalpy_offset;
  const Real liquid_mass = liquid_fraction * densities.liquid;
  const Real solid_mass = (1.0 - liquid_fraction) * densities.solid;
  return {liquid_mass * (h_l + _latent_heat) + solid_mass * h_s,
          liquid_mass * cp_l + solid_mass * cp_s,
          densities.liquid * (h_l + _latent_heat) - densities.solid * h_s};
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::liquidFractionAtTemperature(
    const Real temperature) const
{
  if (_latent_heat <= 0)
    return 1.0;
  if (isothermal())
    return temperature >= _solidus_temperature ? 1.0 : 0.0;
  return liquidFraction(temperature).first;
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::liquidPhaseConductivity(
    const Real temperature, const Real liquid_fraction) const
{
  return liquid_fraction * _fp[LIQUID]->k_from_p_T(_thermodynamic_pressure, temperature) +
         (1.0 - liquid_fraction) * _solid_fp->k_from_p_T(_thermodynamic_pressure, temperature);
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::phaseEnthalpyAtTemperature(
    const Phase phase, const Real temperature) const
{
  if (phase == GAS)
    return phaseEnthalpy(GAS, temperature).first;
  // Per unit mass at equal densities; the kernels only prescribe the enthalpy of inflowing liquid
  const Real liquid_fraction = liquidFractionAtTemperature(temperature);
  return liquidPhaseEnthalpyDensity(temperature, liquid_fraction, {1.0, 1.0, 1.0}).value;
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::phaseDensity(const Phase phase,
                                                          const Moose::FaceArg & face) const
{
  return (*_density[phase])(face, Moose::currentState());
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::enthalpyDensity(const Real alpha,
                                                             const Real temperature,
                                                             const Moose::FaceArg & face) const
{
  const auto state = Moose::currentState();
  const Densities densities{
      phaseDensity(LIQUID, face), (*_solid_density)(face, state), phaseDensity(GAS, face)};
  return alpha * liquidPhaseEnthalpyDensity(
                     temperature, liquidFractionAtTemperature(temperature), densities)
                     .value +
         (1.0 - alpha) * densities.gas * phaseEnthalpyAtTemperature(GAS, temperature);
}

ConservativeDiffuseInterfaceMixtureEnthalpy::CellLinearization
ConservativeDiffuseInterfaceMixtureEnthalpy::cellLinearization(const ElemInfo & elem_info) const
{
  const auto state = Moose::currentState();
  const auto elem_arg = makeElemArg(elem_info.elem());
  const Real alpha = cellAlpha(elem_info);
  const std::array<Real, num_phases> fraction{alpha, 1.0 - alpha};
  const Real liquid_fraction =
      _liquid_fraction_var ? _liquid_fraction_var->getElemValue(elem_info, state) : 0.0;

  CellLinearization lin{_temperature_var->getElemValue(elem_info, state), 0.0, 0.0, {}, {}, 0.0};
  const auto densities = cellDensities(elem_arg, state);
  const Real alpha_rho_gas = fraction[GAS] * densities.gas;

  const auto liquid = liquidPhaseEnthalpyDensity(lin.temperature, liquid_fraction, densities);
  // What the liquid phase advects is its liquid, the solid being held at rest: the liquid that a
  // freezing cell expels as it expands carries the latent heat away with it, at the temperature
  // it left the cell at, which is not below the solidus. A cell that froze within the step and
  // is colder than that then passes the liquid expelled below it on unchanged.
  const bool advects_at_cell_temperature =
      _latent_heat <= 0 || lin.temperature >= _solidus_temperature;
  const auto [h_l, cp_l] =
      phaseEnthalpy(LIQUID, advects_at_cell_temperature ? lin.temperature : _solidus_temperature);
  const auto [h_g, cp_g] = phaseEnthalpy(GAS, lin.temperature);
  lin.phase_enthalpy = {h_l + _latent_heat, h_g};
  lin.enthalpy_density = alpha * liquid.value + alpha_rho_gas * h_g;
  const Real heat_capacity = alpha * liquid.temperature_derivative + alpha_rho_gas * cp_g;

  // How the liquid fraction moves with rho*h: not at all outside the melting interval; with the
  // temperature inside a finite one; almost alone at the melting temperature of an isothermal
  // melt, which the recovery returns exactly, so that the cell stays there while it changes phase
  Real liquid_fraction_derivative = 0.0;
  lin.temperature_derivative = 1.0 / heat_capacity;
  if (hasPhaseChange(alpha * densities.liquid, lin.enthalpy_density))
  {
    const Real latent_capacity = alpha * liquid.liquid_fraction_derivative;
    if (!isothermal())
    {
      const Real df_dT = liquidFraction(lin.temperature).second;
      lin.temperature_derivative = 1.0 / (heat_capacity + latent_capacity * df_dT);
      liquid_fraction_derivative = df_dT * lin.temperature_derivative;
    }
    else if (lin.temperature == _solidus_temperature)
    {
      const Real df_dT = 1.0 / isothermal_linearization_width;
      lin.temperature_derivative = 1.0 / (heat_capacity + latent_capacity * df_dT);
      liquid_fraction_derivative = df_dT * lin.temperature_derivative;
    }
  }
  // In a cell changing phase the advected enthalpy is linearized as (h_l + L) f_l / f_l^*, which
  // has the same value at convergence and lets the outflow fall with the liquid left: a face flux
  // lagged from the previous outer iteration then cannot drain the cell of more latent heat than
  // it holds. A cell below the solidus advects the solidus enthalpy explicitly.
  Real liquid_advected_derivative = 0.0;
  if (advects_at_cell_temperature)
    liquid_advected_derivative =
        (liquid_fraction_derivative > 0 && liquid_fraction > 0)
            ? lin.phase_enthalpy[LIQUID] * liquid_fraction_derivative / liquid_fraction
            : cp_l * lin.temperature_derivative;
  lin.phase_enthalpy_derivative = {liquid_advected_derivative, cp_g * lin.temperature_derivative};

  // In series, the phases conduct across the interface: summed over the cells of a diffuse
  // interface, the cell resistances then add up to the two layers' own, sum_k H_k / k_k, however
  // many cells the interface spreads over, whereas volume-weighting k lets heat leak through it.
  // In parallel, a cell holding a little of the other phase keeps most of its conductivity, which
  // is what a smeared interface along the heat flux needs.
  const Real k_liquid_phase = liquidPhaseConductivity(lin.temperature, liquid_fraction);
  const Real k_gas = _fp[GAS]->k_from_p_T(_thermodynamic_pressure, lin.temperature);
  lin.conductivity = _parallel_conductivity
                         ? fraction[LIQUID] * k_liquid_phase + fraction[GAS] * k_gas
                         : 1.0 / (fraction[LIQUID] / k_liquid_phase + fraction[GAS] / k_gas);
  return lin;
}

ConservativeDiffuseInterfaceMixtureEnthalpy::EnthalpyWithDerivatives
ConservativeDiffuseInterfaceMixtureEnthalpy::cellEnthalpy(const Real alpha,
                                                          const Densities & densities,
                                                          const Real temperature,
                                                          const Real liquid_fraction) const
{
  const auto liquid = liquidPhaseEnthalpyDensity(temperature, liquid_fraction, densities);
  const auto [h_g, cp_g] = phaseEnthalpy(GAS, temperature);
  const Real alpha_rho_gas = (1.0 - alpha) * densities.gas;
  return {alpha * liquid.value + alpha_rho_gas * h_g,
          alpha * liquid.temperature_derivative + alpha_rho_gas * cp_g,
          alpha * liquid.liquid_fraction_derivative};
}

Real
ConservativeDiffuseInterfaceMixtureEnthalpy::solveTemperature(const Real alpha,
                                                              const Densities & densities,
                                                              const Real enthalpy_density,
                                                              const Real initial_temperature,
                                                              const Real T_min,
                                                              const Real T_max,
                                                              const bool mushy,
                                                              const Real fixed_liquid_fraction,
                                                              const Elem & elem) const
{
  // rho*h is monotonic in T on every branch, and the classification puts the root inside
  // [T_min, T_max], so the residual sign brackets it. The slope of a mushy branch can be far from
  // its secant (the linear liquid fraction has no slope at the ends of the interval, the
  // smoothstep little), so a Newton step that leaves the bracket is replaced by bisection.
  Real lower = T_min;
  Real upper = T_max;
  Real temperature = std::clamp(initial_temperature, T_min, T_max);
  for (unsigned int it = 0; it < max_newton_iterations; ++it)
  {
    const auto [f, df_dT] =
        mushy ? liquidFraction(temperature) : std::make_pair(fixed_liquid_fraction, 0.0);
    const auto H = cellEnthalpy(alpha, densities, temperature, f);
    const Real residual = H.value - enthalpy_density;
    (residual > 0 ? upper : lower) = temperature;

    Real new_temperature =
        temperature - residual / (H.temperature_derivative + H.liquid_fraction_derivative * df_dT);
    if (!(new_temperature > lower && new_temperature < upper))
      new_temperature = (std::isfinite(lower) && std::isfinite(upper))
                            ? 0.5 * (lower + upper)
                            : std::clamp(new_temperature, T_min, T_max);

    const Real tolerance = temperature_relative_tolerance * std::abs(new_temperature);
    const bool converged =
        std::abs(new_temperature - temperature) <= tolerance || upper - lower <= tolerance;
    temperature = new_temperature;
    if (converged)
      return temperature;
  }
  mooseError("The temperature recovery from rho*h = ",
             enthalpy_density,
             " did not converge in element ",
             elem.id(),
             " (liquid partial density ",
             alpha * densities.liquid,
             ").");
}

void
ConservativeDiffuseInterfaceMixtureEnthalpy::recoverTemperature()
{
  const auto state = Moose::currentState();
  const auto t_sys = _temperature_var->sys().number();
  const auto t_var = _temperature_var->number();
  auto & temperature_solution = _temperature_var->sys().solution();
  constexpr Real infinity = std::numeric_limits<Real>::infinity();

  for (const auto * elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (elem_info->elem()->processor_id() != processor_id() ||
        !hasBlocks(elem_info->subdomain_id()))
      continue;

    const auto elem_arg = makeElemArg(elem_info->elem());
    const Real alpha = cellAlpha(*elem_info);
    const auto densities = cellDensities(elem_arg, state);
    const Real enthalpy_density = _enthalpy_var->getElemValue(*elem_info, state);
    const Real lagged_temperature = _temperature_var->getElemValue(*elem_info, state);
    Real liquid_fraction =
        _liquid_fraction_var ? _liquid_fraction_var->getElemValue(*elem_info, state) : 0.0;

    Real temperature;
    if (!hasPhaseChange(alpha * densities.liquid, enthalpy_density))
      temperature = solveTemperature(alpha,
                                     densities,
                                     enthalpy_density,
                                     lagged_temperature,
                                     -infinity,
                                     infinity,
                                     false,
                                     liquid_fraction,
                                     *elem_info->elem());
    else
    {
      // Classify rho*h against the fully solid and the fully liquid enthalpy of the melting
      // interval; each of the three branches is smooth in T
      const Real solid_enthalpy = cellEnthalpy(alpha, densities, _solidus_temperature, 0.0).value;
      const Real liquid_enthalpy = cellEnthalpy(alpha, densities, _liquidus_temperature, 1.0).value;
      if (enthalpy_density < solid_enthalpy)
      {
        liquid_fraction = 0.0;
        temperature = solveTemperature(alpha,
                                       densities,
                                       enthalpy_density,
                                       lagged_temperature,
                                       -infinity,
                                       _solidus_temperature,
                                       false,
                                       liquid_fraction,
                                       *elem_info->elem());
      }
      else if (enthalpy_density > liquid_enthalpy)
      {
        liquid_fraction = 1.0;
        temperature = solveTemperature(alpha,
                                       densities,
                                       enthalpy_density,
                                       lagged_temperature,
                                       _liquidus_temperature,
                                       infinity,
                                       false,
                                       liquid_fraction,
                                       *elem_info->elem());
      }
      else if (isothermal())
      {
        temperature = _solidus_temperature;
        liquid_fraction = std::clamp(
            (enthalpy_density - solid_enthalpy) / (liquid_enthalpy - solid_enthalpy), 0.0, 1.0);
      }
      else
      {
        temperature = solveTemperature(alpha,
                                       densities,
                                       enthalpy_density,
                                       lagged_temperature,
                                       _solidus_temperature,
                                       _liquidus_temperature,
                                       true,
                                       0.0,
                                       *elem_info->elem());
        liquid_fraction = liquidFraction(temperature).first;
      }
    }

    temperature_solution.set(elem_info->dofIndices()[t_sys][t_var], temperature);
    if (_liquid_fraction_var)
      _liquid_fraction_var->sys().solution().set(
          elem_info
              ->dofIndices()[_liquid_fraction_var->sys().number()][_liquid_fraction_var->number()],
          liquid_fraction);
  }

  temperature_solution.close();
  _temperature_var->sys().system().update();
  if (_liquid_fraction_var && &_liquid_fraction_var->sys() != &_temperature_var->sys())
  {
    _liquid_fraction_var->sys().solution().close();
    _liquid_fraction_var->sys().system().update();
  }
}

void
ConservativeDiffuseInterfaceMixtureEnthalpy::execute()
{
  // A recovered run restores rho*h from the checkpoint
  if (_app.isRecovering())
    return;

  // The initial liquid fraction follows the initial temperature; an isothermal melt starts liquid
  // at its melting temperature
  if (_liquid_fraction_var)
  {
    const auto state = Moose::currentState();
    auto & liquid_fraction_solution = _liquid_fraction_var->sys().solution();
    for (const auto * elem_info : _fe_problem.mesh().elemInfoVector())
    {
      if (elem_info->elem()->processor_id() != processor_id() ||
          !hasBlocks(elem_info->subdomain_id()))
        continue;

      const Real temperature = _temperature_var->getElemValue(*elem_info, state);
      liquid_fraction_solution.set(elem_info->dofIndices()[_liquid_fraction_var->sys().number()]
                                                          [_liquid_fraction_var->number()],
                                   isothermal() ? (temperature >= _solidus_temperature ? 1.0 : 0.0)
                                                : liquidFraction(temperature).first);
    }
    liquid_fraction_solution.close();
    _liquid_fraction_var->sys().system().update();
  }

  const auto h_sys = _enthalpy_var->sys().number();
  const auto h_var = _enthalpy_var->number();
  auto & enthalpy_solution = _enthalpy_var->sys().solution();

  for (const auto * elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (elem_info->elem()->processor_id() != processor_id() ||
        !hasBlocks(elem_info->subdomain_id()))
      continue;

    enthalpy_solution.set(elem_info->dofIndices()[h_sys][h_var],
                          cellLinearization(*elem_info).enthalpy_density);
  }

  enthalpy_solution.close();
  _enthalpy_var->sys().system().update();
}
