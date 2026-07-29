//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

#include <array>

class LaserBeamFoamLiquidFractionCorrector;

/**
 * Declares the low-Mach phase-change divergence source implied by unequal
 * liquid and solid densities:
 *
 *   div(u)_pc = alpha (rho_s - rho_l) / rho * D(fl)/Dt
 */
class EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();
  EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial(const InputParameters & params);

protected:
  const Moose::Functor<ADReal> & _liquid_fraction;
  const Moose::Functor<ADReal> & _density;
  const Moose::Functor<ADReal> & _rho_liquid;
  const Moose::Functor<ADReal> & _rho_solid;
  const Moose::Functor<ADReal> & _material_fraction;
  std::array<const Moose::Functor<ADReal> *, LIBMESH_DIM> _velocity;
  const LaserBeamFoamLiquidFractionCorrector * _liquid_fraction_corrector;
  const Real _min_density;

  template <typename SpaceArg>
  Real liquidFractionRate(const SpaceArg & r, const Moose::StateArg & t) const;

  using UserObjectInterface::getUserObject;
};
