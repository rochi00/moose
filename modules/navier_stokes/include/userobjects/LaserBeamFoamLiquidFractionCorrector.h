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
#include "MooseVariableFV.h"
#include "MooseFunctorArguments.h"

#include <unordered_map>

class SystemBase;

/**
 * Explicit enthalpy liquid-fraction correction.
 *
 * The correction is
 *
 *   f_l <- clamp(f_l + relaxation * cp / L * (T - T_corr), 0, 1)
 *   T_corr = (T_liquidus - T_solidus) * f_l + T_solidus
 *
 * This object intentionally only updates the liquid-fraction field. The
 * temperature equation is still solved by the existing PIMPLE energy system.
 */
class LaserBeamFoamLiquidFractionCorrector : public GeneralUserObject,
                                             public NonADFunctorInterface,
                                             public BlockRestrictable
{
public:
  static InputParameters validParams();

  LaserBeamFoamLiquidFractionCorrector(const InputParameters & params);

  void initialize() override {}
  void execute() override;
  void finalize() override {}

  Real correctLiquidFraction();

  Real maxCorrection() const { return _max_correction; }
  Real meanCorrection() const { return _mean_correction; }
  Real phaseChangeRate(const Elem * elem) const;
  template <typename SpaceArg>
  Real phaseChangeRate(const SpaceArg & arg) const
  {
    return phaseChangeRate(arg.elem);
  }
  Real phaseChangeRate(const Moose::ElemArg & elem_arg) const;
  Real phaseChangeRate(const Moose::FaceArg & face_arg) const;
  Real phaseChangeRate(const Moose::NodeArg &) const { return 0.0; }

protected:
  MooseVariableFVReal & _liquid_fraction_var;
  SystemBase & _liquid_fraction_system;
  const unsigned int _liquid_fraction_sys_num;
  const unsigned int _liquid_fraction_var_num;

  const Moose::Functor<Real> & _temperature;
  const Moose::Functor<Real> & _specific_heat;
  const Moose::Functor<Real> & _latent_heat;
  const Moose::Functor<Real> & _T_solidus;
  const Moose::Functor<Real> & _T_liquidus;

  const Moose::Functor<Real> * const _material_fraction;
  const Real _min_material_fraction;
  const Real _relaxation;
  const Real _max_liquid_fraction_change;
  const Real _phase_change_rate_relaxation;

  std::unordered_map<dof_id_type, Real> _phase_change_rate;

  Real _max_correction = 0.0;
  Real _mean_correction = 0.0;
};
