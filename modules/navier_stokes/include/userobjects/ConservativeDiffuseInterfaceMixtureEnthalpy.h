//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "NonADFunctorInterface.h"
#include "BlockRestrictable.h"
#include "MooseLinearVariableFV.h"

#include <array>

class SinglePhaseFluidProperties;
class ElemInfo;

/**
 * Relates the conserved mixture enthalpy density rho*h = sum_k alpha_k rho_k h_k(T) of a
 * two-phase diffuse-interface flow to its temperature. The phase densities are the constant ones
 * of the volume-fraction equation; the phase enthalpies come from fluid properties evaluated at a
 * fixed thermodynamic pressure. The liquid phase may melt and solidify: its enthalpy then carries
 * the latent heat L f_l of its liquid fraction f_l, which is recovered from rho*h together with the
 * temperature, over a finite melting interval or at a single melting temperature, and its solid
 * may have its own density, heat capacity and conductivity. Provides the
 * per-cell linearization of rho*h about the lagged state that the energy kernels assemble with,
 * recovers the temperature from the solved rho*h, and sets the initial rho*h from the initial
 * temperature and volume fraction.
 */
class ConservativeDiffuseInterfaceMixtureEnthalpy : public GeneralUserObject,
                                                    public NonADFunctorInterface,
                                                    public BlockRestrictable
{
public:
  static InputParameters validParams();

  ConservativeDiffuseInterfaceMixtureEnthalpy(const InputParameters & params);

  void initialSetup() override;
  void meshChanged() override;
  void initialize() override {}
  /// Sets rho*h from the initial temperature and volume fraction
  void execute() override;
  void finalize() override {}

  /// Indices of the two phases; the volume fraction is the liquid's
  enum Phase : unsigned int
  {
    LIQUID = 0,
    GAS = 1
  };
  static constexpr unsigned int num_phases = 2;

  /// Linearization of rho*h about the lagged temperature and liquid fraction of a cell, at the
  /// current volume fraction
  struct CellLinearization
  {
    /// The lagged temperature
    Real temperature;
    /// rho*h of the current volume fraction at the lagged state
    Real enthalpy_density;
    /// dT/d(rho*h) at the lagged state: 1 / sum_k alpha_k rho_k cp_k outside the melting interval,
    /// smaller inside it, and nearly zero at the melting temperature of an isothermal melt
    Real temperature_derivative;
    /// Specific enthalpy h_k that the mass flux of each phase advects at the lagged state; for the
    /// liquid phase that of its liquid, h_l + L, its solid being at rest
    std::array<Real, num_phases> phase_enthalpy;
    /// dh_k/d(rho*h) of each phase at the lagged state
    std::array<Real, num_phases> phase_enthalpy_derivative;
    /// Mixture thermal conductivity at the lagged temperature, the series mix 1/k = sum_k alpha_k / k_k
    Real conductivity;
  };
  CellLinearization cellLinearization(const ElemInfo & elem_info) const;

  /// An enthalpy with its derivatives in temperature and in the liquid fraction
  struct EnthalpyWithDerivatives
  {
    Real value;
    Real temperature_derivative;
    Real liquid_fraction_derivative;
  };

  /// The lagged temperature of a cell expressed through rho*h:
  /// T = rho*h dT/d(rho*h) + (T_lag - (rho*h)^* dT/d(rho*h)), returning the constant in parentheses
  static Real temperatureOffset(const CellLinearization & lin)
  {
    return lin.temperature - lin.enthalpy_density * lin.temperature_derivative;
  }

  /// The densities of the liquid phase's liquid and solid, and of the gas
  struct Densities
  {
    Real liquid;
    Real solid;
    Real gas;
  };

  /// Sensible specific enthalpy of a phase's fluid properties and its temperature derivative cp
  std::pair<Real, Real> phaseEnthalpy(Phase phase, Real temperature) const;

  /// Enthalpy per unit volume of the liquid phase at a temperature and liquid fraction f, the
  /// volume fraction of its liquid: f rho_l [h_l(T) + L] + (1 - f) rho_s h_s(T), with the solid
  /// enthalpy h_s referenced to the liquid's at the solidus so that the two differ by the latent
  /// heat alone there
  EnthalpyWithDerivatives liquidPhaseEnthalpyDensity(Real temperature,
                                                     Real liquid_fraction,
                                                     const Densities & densities) const;

  /// Specific enthalpy of a phase at a temperature, including the latent heat of the liquid fraction
  /// that temperature implies (fully liquid at the melting temperature of an isothermal melt)
  Real phaseEnthalpyAtTemperature(Phase phase, Real temperature) const;

  /// rho*h of a volume fraction at a temperature
  Real enthalpyDensity(Real alpha, Real temperature, const Moose::FaceArg & face) const;

  /// The constant density of a phase on a face
  Real phaseDensity(Phase phase, const Moose::FaceArg & face) const;

  /// Overwrites the temperature, and the liquid fraction when there is latent heat, with the root
  /// of alpha {f_l rho_l [h_l(T) + L] + (1 - f_l) rho_s h_s(T)} + (1 - alpha) rho_g h_g(T) = rho*h
  /// in every cell
  void recoverTemperature();

  const VariableName & enthalpyVariableName() const { return _enthalpy_name; }

private:
  void cacheVariables();
  Real cellAlpha(const ElemInfo & elem_info) const;
  Densities cellDensities(const Moose::ElemArg & elem_arg, const Moose::StateArg & state) const;
  /// The liquid fraction a temperature implies, fully liquid at the melting temperature of an
  /// isothermal melt, and immaterial without latent heat
  Real liquidFractionAtTemperature(Real temperature) const;
  /// Whether the melting interval has zero width
  bool isothermal() const { return _liquidus_temperature == _solidus_temperature; }
  /// Liquid fraction of a finite melting interval at a temperature, and its temperature derivative
  std::pair<Real, Real> liquidFraction(Real temperature) const;
  /// Whether a cell of liquid partial density alpha rho_l and enthalpy scale H changes phase: its
  /// latent heat alpha rho_l L must be resolvable against the round-off of rho*h
  bool hasPhaseChange(Real liquid_partial_density, Real enthalpy_scale) const;
  /// Thermal conductivity of the liquid phase at a temperature and liquid fraction: its liquid and
  /// solid conduct in parallel, f k_l + (1 - f) k_s
  Real liquidPhaseConductivity(Real temperature, Real liquid_fraction) const;
  /// rho*h of a cell at a temperature and liquid fraction, with its derivatives in each
  EnthalpyWithDerivatives cellEnthalpy(Real alpha,
                                       const Densities & densities,
                                       Real temperature,
                                       Real liquid_fraction) const;
  /// Newton solve of rho*h = H(T) on one branch of the enthalpy, T clamped to [T_min, T_max]. The
  /// liquid fraction is held at fixed_liquid_fraction, or follows liquidFraction(T) on the mushy
  /// branch of a finite melting interval.
  Real solveTemperature(Real alpha,
                        const Densities & densities,
                        Real enthalpy_density,
                        Real initial_temperature,
                        Real T_min,
                        Real T_max,
                        bool mushy,
                        Real fixed_liquid_fraction,
                        const Elem & elem) const;

  const VariableName _enthalpy_name;
  const VariableName _temperature_name;
  const Moose::Functor<Real> & _alpha;
  const std::array<const Moose::Functor<Real> *, num_phases> _density;
  /// Density of the liquid phase's solid; the liquid's own without solid_density
  const Moose::Functor<Real> * const _solid_density;
  std::array<const SinglePhaseFluidProperties *, num_phases> _fp;
  /// Fluid properties of the liquid phase's solid; the liquid's own without solid_fluid_properties
  const SinglePhaseFluidProperties * const _solid_fp;
  /// Added to the solid enthalpy so that it equals the liquid's at the solidus
  Real _solid_enthalpy_offset = 0.0;
  /// Pressure at which the phase enthalpies are evaluated. The flow solves a reduced pressure,
  /// which is not a thermodynamic pressure, and the enthalpy of a low-Mach flow does not depend on
  /// its pressure variations.
  const Real _thermodynamic_pressure;

  /// Latent heat of fusion of the liquid phase [J/kg]; zero without phase change
  const Real _latent_heat;
  const Real _solidus_temperature;
  const Real _liquidus_temperature;
  /// Whether the liquid fraction is a smoothstep rather than linear over the melting interval
  const bool _smoothstep;
  /// Whether the phase conductivities combine in parallel rather than in series in a cell
  const bool _parallel_conductivity;

  MooseLinearVariableFVReal * _enthalpy_var = nullptr;
  MooseLinearVariableFVReal * _temperature_var = nullptr;
  /// The liquid fraction of the liquid phase, recovered from rho*h; null without latent heat
  MooseLinearVariableFVReal * _liquid_fraction_var = nullptr;
};
