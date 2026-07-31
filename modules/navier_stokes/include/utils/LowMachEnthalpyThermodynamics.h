//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

#include <algorithm>
#include <cmath>

namespace LowMachEnthalpy
{
struct Parameters
{
  Real cp_gas;
  Real cp_solid;
  Real cp_liquid;
  Real rho_gas;
  Real rho_solid;
  Real rho_liquid;
  Real latent_heat;
  Real T_solidus;
  Real T_liquidus;
  Real T_reference;
};

inline Real
solidEnthalpy(const Parameters & p)
{
  return p.cp_solid * (p.T_solidus - p.T_reference);
}

inline Real
liquidEnthalpy(const Parameters & p)
{
  return 0.5 * (p.cp_solid + p.cp_liquid) * (p.T_liquidus - p.T_solidus) + solidEnthalpy(p) +
         p.latent_heat;
}

inline bool
isPCM(const Real material_fraction)
{
  return material_fraction >= 0.5;
}

inline Real
enthalpyFromTemperature(const Real temperature, const Real material_fraction, const Parameters & p)
{
  if (!isPCM(material_fraction))
    return p.cp_gas * (temperature - p.T_reference);

  const Real h_solid = solidEnthalpy(p);
  const Real h_liquid = liquidEnthalpy(p);

  if (temperature < p.T_solidus)
    return p.cp_solid * (temperature - p.T_reference);
  if (temperature <= p.T_liquidus)
    return h_solid +
           (h_liquid - h_solid) * (temperature - p.T_solidus) / (p.T_liquidus - p.T_solidus);
  return h_liquid + p.cp_liquid * (temperature - p.T_liquidus);
}

inline Real
temperatureFromEnthalpy(const Real enthalpy, const Real material_fraction, const Parameters & p)
{
  if (!isPCM(material_fraction))
    return enthalpy / p.cp_gas + p.T_reference;

  const Real h_solid = solidEnthalpy(p);
  const Real h_liquid = liquidEnthalpy(p);

  if (enthalpy < h_solid)
    return enthalpy / p.cp_solid + p.T_reference;
  if (enthalpy <= h_liquid)
    return p.T_solidus + (enthalpy - h_solid) * (p.T_liquidus - p.T_solidus) / (h_liquid - h_solid);
  return p.T_liquidus + (enthalpy - h_liquid) / p.cp_liquid;
}

inline Real
liquidFractionFromEnthalpy(const Real enthalpy, const Real material_fraction, const Parameters & p)
{
  if (!isPCM(material_fraction))
    return 0.0;

  const Real h_solid = solidEnthalpy(p);
  const Real h_liquid = liquidEnthalpy(p);

  if (enthalpy < h_solid)
    return 0.0;
  if (enthalpy > h_liquid)
    return 1.0;

  const Real denominator =
      enthalpy * (p.rho_liquid - p.rho_solid) - p.rho_liquid * h_liquid + p.rho_solid * h_solid;
  return std::clamp(p.rho_solid * (h_solid - enthalpy) / denominator, 0.0, 1.0);
}

inline Real
liquidFractionFromEnthalpy(const Real enthalpy,
                           const Real material_fraction,
                           const Real rho_solid,
                           const Real rho_liquid,
                           const Real h_solid,
                           const Real h_liquid)
{
  if (!isPCM(material_fraction))
    return 0.0;

  if (enthalpy < h_solid)
    return 0.0;
  if (enthalpy > h_liquid)
    return 1.0;

  const Real denominator =
      enthalpy * (rho_liquid - rho_solid) - rho_liquid * h_liquid + rho_solid * h_solid;
  return std::clamp(rho_solid * (h_solid - enthalpy) / denominator, 0.0, 1.0);
}

inline Real
dEnthalpyDTemperature(const Real temperature, const Real material_fraction, const Parameters & p)
{
  if (!isPCM(material_fraction))
    return p.cp_gas;
  if (temperature < p.T_solidus)
    return p.cp_solid;
  if (temperature <= p.T_liquidus)
    return 0.5 * (p.cp_solid + p.cp_liquid) + p.latent_heat / (p.T_liquidus - p.T_solidus);
  return p.cp_liquid;
}

inline Real
newtonUpdatedEnthalpy(const Real old_enthalpy,
                      const Real old_dh_dT,
                      const Real old_temperature,
                      const Real linearized_temperature)
{
  return old_enthalpy + old_dh_dT * (linearized_temperature - old_temperature);
}

inline Real
dLiquidFractionDEnthalpy(const Real enthalpy, const Real material_fraction, const Parameters & p)
{
  if (!isPCM(material_fraction))
    return 0.0;

  const Real h_solid = solidEnthalpy(p);
  const Real h_liquid = liquidEnthalpy(p);
  if (enthalpy < h_solid || enthalpy > h_liquid)
    return 0.0;

  const Real denominator =
      enthalpy * (p.rho_liquid - p.rho_solid) - p.rho_liquid * h_liquid + p.rho_solid * h_solid;
  return -p.rho_solid * p.rho_liquid * (h_solid - h_liquid) / (denominator * denominator);
}

inline Real
density(const Real material_fraction, const Real liquid_fraction, const Parameters & p)
{
  return p.rho_gas + (p.rho_solid - p.rho_gas) * material_fraction +
         (p.rho_liquid - p.rho_solid) * material_fraction * liquid_fraction;
}

inline Real
threePhaseProperty(const Real material_fraction,
                   const Real liquid_fraction,
                   const Real gas_value,
                   const Real solid_value,
                   const Real liquid_value)
{
  return gas_value + (solid_value - gas_value) * material_fraction +
         (liquid_value - solid_value) * material_fraction * liquid_fraction;
}

inline Real
divergenceSource(const Real material_fraction,
                 const Real eos_density,
                 const Real rho_solid,
                 const Real rho_liquid,
                 const Real dliquid_fraction_dh,
                 const Real heating_rate)
{
  return (rho_solid - rho_liquid) * material_fraction * dliquid_fraction_dh * heating_rate /
         (eos_density * eos_density);
}

/**
 * Backward-Euler phase-change divergence compatible with implicit conservative mass transport.
 *
 * The conservative enthalpy balance supplies the end-of-step volumetric heating rate q. Backward
 * Euler gives h_old_material = h_new - dt q / rho_new. The exact enthalpy-to-liquid-fraction map
 * then gives the density before the thermodynamic update. For
 *
 *   (rho_new - rho_old_material) / dt + rho_new div(u_new) = 0,
 *
 * choosing div(u_new) = (rho_old_material / rho_new - 1) / dt makes the implicit mass update
 * exactly compatible with the finite EOS density change.
 */
inline Real
backwardEulerDivergenceSource(const Real material_fraction,
                              const Real eos_density,
                              const Real rho_solid,
                              const Real rho_liquid,
                              const Real current_enthalpy,
                              const Real h_solid,
                              const Real h_liquid,
                              const Real heating_rate,
                              const Real dt)
{
  if (!isPCM(material_fraction) || rho_solid == rho_liquid || heating_rate == 0.0)
    return 0.0;

  const Real previous_enthalpy = current_enthalpy - dt * heating_rate / eos_density;
  const Real current_liquid_fraction = liquidFractionFromEnthalpy(
      current_enthalpy, material_fraction, rho_solid, rho_liquid, h_solid, h_liquid);
  const Real previous_liquid_fraction = liquidFractionFromEnthalpy(
      previous_enthalpy, material_fraction, rho_solid, rho_liquid, h_solid, h_liquid);
  const Real previous_density =
      eos_density + (rho_liquid - rho_solid) * material_fraction *
                        (previous_liquid_fraction - current_liquid_fraction);

  return (previous_density / eos_density - 1.0) / dt;
}

inline Real
solidDragCoefficient(const Real material_fraction,
                     const Real liquid_fraction,
                     const Real rho_solid,
                     const Real dt,
                     const Real regularization)
{
  const Real solid_fraction = material_fraction * (1.0 - liquid_fraction);
  return rho_solid / dt * solid_fraction * solid_fraction /
         (std::pow(1.0 - solid_fraction, 3) + regularization);
}

inline Real
volumetricLatentEnthalpy(const Real material_fraction,
                         const Real liquid_fraction,
                         const Parameters & p)
{
  return material_fraction * p.rho_liquid * p.latent_heat * liquid_fraction;
}
}
