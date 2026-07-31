//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LowMachEnthalpyThermodynamicsFunctorMaterial.h"

#include "LowMachEnthalpyThermodynamics.h"

registerMooseObject("NavierStokesApp", LowMachEnthalpyThermodynamicsFunctorMaterial);

InputParameters
LowMachEnthalpyThermodynamicsFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Defines the exact enthalpy-temperature, liquid-fraction, and unequal-density EOS closure "
      "for the low-Mach gas/liquid/solid phase-change method.");

  params.addRequiredParam<MooseFunctorName>("material_fraction", "Gas/PCM material fraction.");
  params.addRequiredParam<MooseFunctorName>("temperature", "Temperature.");
  params.addRequiredParam<MooseFunctorName>("cp_gas", "Gas specific heat.");
  params.addRequiredParam<MooseFunctorName>("cp_solid", "Solid PCM specific heat.");
  params.addRequiredParam<MooseFunctorName>("cp_liquid", "Liquid PCM specific heat.");
  params.addRequiredParam<MooseFunctorName>("rho_gas", "Gas density.");
  params.addRequiredParam<MooseFunctorName>("rho_solid", "Solid PCM density.");
  params.addRequiredParam<MooseFunctorName>("rho_liquid", "Liquid PCM density.");
  params.addRequiredParam<MooseFunctorName>("k_gas", "Gas thermal conductivity.");
  params.addRequiredParam<MooseFunctorName>("k_solid", "Solid PCM thermal conductivity.");
  params.addRequiredParam<MooseFunctorName>("k_liquid", "Liquid PCM thermal conductivity.");
  params.addRequiredParam<MooseFunctorName>("mu_gas", "Gas dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("mu_solid", "Solid PCM dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("mu_liquid", "Liquid PCM dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("L", "Latent heat.");
  params.addRequiredParam<MooseFunctorName>("T_solidus", "Solidus temperature.");
  params.addRequiredParam<MooseFunctorName>("T_liquidus", "Liquidus temperature.");
  params.addParam<MooseFunctorName>("T_reference", 0.0, "Reference temperature.");

  params.addParam<std::string>("enthalpy_from_temperature_name",
                               "low_mach_enthalpy_from_temperature",
                               "Name of the specific-enthalpy-from-temperature functor.");
  params.addParam<std::string>("temperature_from_enthalpy_name",
                               "low_mach_temperature_from_enthalpy",
                               "Name of the temperature-from-specific-enthalpy functor.");
  params.addParam<std::string>(
      "liquid_fraction_name", "low_mach_liquid_fraction", "Name of the liquid-fraction functor.");
  params.addParam<std::string>(
      "density_name", "low_mach_density", "Name of the three-phase EOS density functor.");
  params.addParam<std::string>("specific_heat_name",
                               "low_mach_specific_heat",
                               "Name of the three-phase specific-heat functor.");
  params.addParam<std::string>("thermal_conductivity_name",
                               "low_mach_thermal_conductivity",
                               "Name of the three-phase thermal-conductivity functor.");
  params.addParam<std::string>("dynamic_viscosity_name",
                               "low_mach_dynamic_viscosity",
                               "Name of the three-phase dynamic-viscosity functor.");
  params.addParam<std::string>("solid_drag_name",
                               "low_mach_solid_drag",
                               "Name of the Carman-Kozeny solid-drag coefficient functor.");
  params.addRangeCheckedParam<Real>("solid_drag_regularization",
                                    1e-3,
                                    "solid_drag_regularization>0",
                                    "Carman-Kozeny denominator regularization.");
  params.addParam<std::string>(
      "dh_dT_name", "low_mach_dh_dT", "Name of the enthalpy temperature-derivative functor.");
  params.addParam<std::string>("dliquid_fraction_dh_name",
                               "low_mach_dliquid_fraction_dh",
                               "Name of the liquid-fraction enthalpy-derivative functor.");
  params.addParam<std::string>("volumetric_latent_enthalpy_name",
                               "low_mach_volumetric_latent_enthalpy",
                               "Name of the volumetric latent-enthalpy functor.");
  params.addParam<std::string>(
      "solid_enthalpy_name", "low_mach_solid_enthalpy", "Name of the solidus enthalpy functor.");
  params.addParam<std::string>(
      "liquid_enthalpy_name", "low_mach_liquid_enthalpy", "Name of the liquidus enthalpy functor.");
  return params;
}

LowMachEnthalpyThermodynamicsFunctorMaterial::LowMachEnthalpyThermodynamicsFunctorMaterial(
    const InputParameters & params)
  : FunctorMaterial(params),
    _material_fraction(getFunctor<Real>("material_fraction")),
    _temperature(getFunctor<Real>("temperature")),
    _cp_gas(getFunctor<Real>("cp_gas")),
    _cp_solid(getFunctor<Real>("cp_solid")),
    _cp_liquid(getFunctor<Real>("cp_liquid")),
    _rho_gas(getFunctor<Real>("rho_gas")),
    _rho_solid(getFunctor<Real>("rho_solid")),
    _rho_liquid(getFunctor<Real>("rho_liquid")),
    _k_gas(getFunctor<Real>("k_gas")),
    _k_solid(getFunctor<Real>("k_solid")),
    _k_liquid(getFunctor<Real>("k_liquid")),
    _mu_gas(getFunctor<Real>("mu_gas")),
    _mu_solid(getFunctor<Real>("mu_solid")),
    _mu_liquid(getFunctor<Real>("mu_liquid")),
    _latent_heat(getFunctor<Real>("L")),
    _T_solidus(getFunctor<Real>("T_solidus")),
    _T_liquidus(getFunctor<Real>("T_liquidus")),
    _T_reference(getFunctor<Real>("T_reference"))
{
  const auto thermodynamic_parameters = [this](const auto & r,
                                               const auto & state) -> LowMachEnthalpy::Parameters
  {
    LowMachEnthalpy::Parameters p{_cp_gas(r, state),
                                  _cp_solid(r, state),
                                  _cp_liquid(r, state),
                                  _rho_gas(r, state),
                                  _rho_solid(r, state),
                                  _rho_liquid(r, state),
                                  _latent_heat(r, state),
                                  _T_solidus(r, state),
                                  _T_liquidus(r, state),
                                  _T_reference(r, state)};
    if (p.cp_gas <= 0.0 || p.cp_solid <= 0.0 || p.cp_liquid <= 0.0)
      mooseError(name(), ": all specific heats must be positive.");
    if (p.rho_gas <= 0.0 || p.rho_solid <= 0.0 || p.rho_liquid <= 0.0)
      mooseError(name(), ": all densities must be positive.");
    if (p.latent_heat < 0.0)
      mooseError(name(), ": latent heat must be nonnegative.");
    if (p.T_liquidus <= p.T_solidus)
      mooseError(name(), ": T_liquidus must be greater than T_solidus.");
    return p;
  };

  addFunctorProperty<Real>(
      getParam<std::string>("enthalpy_from_temperature_name"),
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
      {
        return LowMachEnthalpy::enthalpyFromTemperature(_temperature(r, state),
                                                        _material_fraction(r, state),
                                                        thermodynamic_parameters(r, state));
      });

  const auto liquid_fraction_from_temperature =
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
  {
    const auto p = thermodynamic_parameters(r, state);
    const Real material_fraction = _material_fraction(r, state);
    const Real enthalpy =
        LowMachEnthalpy::enthalpyFromTemperature(_temperature(r, state), material_fraction, p);
    return LowMachEnthalpy::liquidFractionFromEnthalpy(enthalpy, material_fraction, p);
  };

  addFunctorProperty<Real>(
      getParam<std::string>("temperature_from_enthalpy_name"),
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
      {
        const auto p = thermodynamic_parameters(r, state);
        const Real material_fraction = _material_fraction(r, state);
        const Real enthalpy =
            LowMachEnthalpy::enthalpyFromTemperature(_temperature(r, state), material_fraction, p);
        return LowMachEnthalpy::temperatureFromEnthalpy(enthalpy, material_fraction, p);
      });

  addFunctorProperty<Real>(getParam<std::string>("liquid_fraction_name"),
                           liquid_fraction_from_temperature);

  addFunctorProperty<Real>(
      getParam<std::string>("density_name"),
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
      {
        const auto p = thermodynamic_parameters(r, state);
        const Real material_fraction = _material_fraction(r, state);
        const Real enthalpy =
            LowMachEnthalpy::enthalpyFromTemperature(_temperature(r, state), material_fraction, p);
        const Real liquid_fraction =
            LowMachEnthalpy::liquidFractionFromEnthalpy(enthalpy, material_fraction, p);
        return LowMachEnthalpy::density(material_fraction, liquid_fraction, p);
      });

  addFunctorProperty<Real>(
      getParam<std::string>("specific_heat_name"),
      [this, liquid_fraction_from_temperature](const auto & r, const auto & state) -> Real
      {
        return LowMachEnthalpy::threePhaseProperty(_material_fraction(r, state),
                                                   liquid_fraction_from_temperature(r, state),
                                                   _cp_gas(r, state),
                                                   _cp_solid(r, state),
                                                   _cp_liquid(r, state));
      });

  addFunctorProperty<Real>(
      getParam<std::string>("thermal_conductivity_name"),
      [this, liquid_fraction_from_temperature](const auto & r, const auto & state) -> Real
      {
        return LowMachEnthalpy::threePhaseProperty(_material_fraction(r, state),
                                                   liquid_fraction_from_temperature(r, state),
                                                   _k_gas(r, state),
                                                   _k_solid(r, state),
                                                   _k_liquid(r, state));
      });

  addFunctorProperty<Real>(
      getParam<std::string>("dynamic_viscosity_name"),
      [this, liquid_fraction_from_temperature](const auto & r, const auto & state) -> Real
      {
        return LowMachEnthalpy::threePhaseProperty(_material_fraction(r, state),
                                                   liquid_fraction_from_temperature(r, state),
                                                   _mu_gas(r, state),
                                                   _mu_solid(r, state),
                                                   _mu_liquid(r, state));
      });

  const Real solid_drag_regularization = getParam<Real>("solid_drag_regularization");
  addFunctorProperty<Real>(
      getParam<std::string>("solid_drag_name"),
      [this, liquid_fraction_from_temperature, solid_drag_regularization](
          const auto & r, const auto & state) -> Real
      {
        if (_dt <= 0.0)
          mooseError(name(), ": a positive time step is required for the solid drag coefficient.");
        return LowMachEnthalpy::solidDragCoefficient(_material_fraction(r, state),
                                                     liquid_fraction_from_temperature(r, state),
                                                     _rho_solid(r, state),
                                                     _dt,
                                                     solid_drag_regularization);
      });

  addFunctorProperty<Real>(
      getParam<std::string>("dh_dT_name"),
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
      {
        return LowMachEnthalpy::dEnthalpyDTemperature(_temperature(r, state),
                                                      _material_fraction(r, state),
                                                      thermodynamic_parameters(r, state));
      });

  addFunctorProperty<Real>(
      getParam<std::string>("dliquid_fraction_dh_name"),
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
      {
        const auto p = thermodynamic_parameters(r, state);
        const Real material_fraction = _material_fraction(r, state);
        const Real enthalpy =
            LowMachEnthalpy::enthalpyFromTemperature(_temperature(r, state), material_fraction, p);
        return LowMachEnthalpy::dLiquidFractionDEnthalpy(enthalpy, material_fraction, p);
      });

  addFunctorProperty<Real>(
      getParam<std::string>("volumetric_latent_enthalpy_name"),
      [this, thermodynamic_parameters](const auto & r, const auto & state) -> Real
      {
        const auto p = thermodynamic_parameters(r, state);
        const Real material_fraction = _material_fraction(r, state);
        const Real enthalpy =
            LowMachEnthalpy::enthalpyFromTemperature(_temperature(r, state), material_fraction, p);
        const Real liquid_fraction =
            LowMachEnthalpy::liquidFractionFromEnthalpy(enthalpy, material_fraction, p);
        return LowMachEnthalpy::volumetricLatentEnthalpy(material_fraction, liquid_fraction, p);
      });

  addFunctorProperty<Real>(
      getParam<std::string>("solid_enthalpy_name"),
      [thermodynamic_parameters](const auto & r, const auto & state) -> Real
      { return LowMachEnthalpy::solidEnthalpy(thermodynamic_parameters(r, state)); });

  addFunctorProperty<Real>(
      getParam<std::string>("liquid_enthalpy_name"),
      [thermodynamic_parameters](const auto & r, const auto & state) -> Real
      { return LowMachEnthalpy::liquidEnthalpy(thermodynamic_parameters(r, state)); });
}
