//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial.h"
#include "LaserBeamFoamLiquidFractionCorrector.h"

#include <cmath>

registerMooseObject("NavierStokesApp", EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial);

InputParameters
EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Declares a phase-change velocity-divergence functor for unequal liquid and solid "
      "densities: material_fraction * (rho_solid - rho_liquid) / rho * D(liquid_fraction)/Dt.");
  params.addParam<std::string>(
      "property_name", "phase_change_divergence", "Name of the declared divergence functor.");
  params.addRequiredParam<MooseFunctorName>("liquid_fraction", "Liquid-fraction functor.");
  params.addRequiredParam<MooseFunctorName>("density", "Mixture density functor.");
  params.addRequiredParam<MooseFunctorName>("rho_liquid", "Liquid density functor.");
  params.addRequiredParam<MooseFunctorName>("rho_solid", "Solid density functor.");
  params.addParam<MooseFunctorName>(
      "material_fraction", 1.0, "Material volume fraction multiplying the phase-change source.");
  params.addParam<std::vector<MooseFunctorName>>(
      "velocity",
      {},
      "Optional velocity-component functors. When supplied, D(liquid_fraction)/Dt includes "
      "u dot grad(liquid_fraction); otherwise only the local time derivative is used.");
  params.addParam<UserObjectName>(
      "liquid_fraction_corrector",
      "Optional liquid-fraction corrector supplying a stored, limited local phase-change rate. "
      "When supplied, this is used instead of differentiating the liquid-fraction functor.");
  params.addRangeCheckedParam<Real>(
      "min_density", 1e-12, "min_density>0", "Minimum absolute density used in the denominator.");
  return params;
}

EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial::
    EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial(const InputParameters & params)
  : FunctorMaterial(params),
    _liquid_fraction(getFunctor<ADReal>("liquid_fraction")),
    _density(getFunctor<ADReal>("density")),
    _rho_liquid(getFunctor<ADReal>("rho_liquid")),
    _rho_solid(getFunctor<ADReal>("rho_solid")),
    _material_fraction(getFunctor<ADReal>("material_fraction")),
    _velocity(),
    _liquid_fraction_corrector(
        isParamValid("liquid_fraction_corrector")
            ? &getUserObject<LaserBeamFoamLiquidFractionCorrector>("liquid_fraction_corrector")
            : nullptr),
    _min_density(getParam<Real>("min_density"))
{
  const auto & velocity_names = getParam<std::vector<MooseFunctorName>>("velocity");
  if (velocity_names.size() > LIBMESH_DIM)
    paramError("velocity", "At most ", LIBMESH_DIM, " velocity components may be supplied.");

  for (const auto i : index_range(velocity_names))
    _velocity[i] = &getFunctor<ADReal>(velocity_names[i]);

  const auto property_name = getParam<std::string>("property_name");
  addFunctorProperty<Real>(property_name,
                           [this](const auto & r, const auto & t) -> Real
                           {
                             const Real rho = MetaPhysicL::raw_value(_density(r, t));
                             const Real abs_rho = std::abs(rho);
                             const Real denominator =
                                 abs_rho > _min_density ? rho : std::copysign(_min_density, rho);

                             const Real liquid_fraction_material_derivative =
                                 liquidFractionRate(r, t);

                             return MetaPhysicL::raw_value(_material_fraction(r, t) *
                                                           (_rho_solid(r, t) - _rho_liquid(r, t))) /
                                    denominator * liquid_fraction_material_derivative;
                           });
}

template <typename SpaceArg>
Real
EnthalpyLiquidFractionPhaseChangeDivergenceFunctorMaterial::liquidFractionRate(
    const SpaceArg & r, const Moose::StateArg & t) const
{
  if (_liquid_fraction_corrector)
    return _liquid_fraction_corrector->phaseChangeRate(r);

  Real liquid_fraction_material_derivative = MetaPhysicL::raw_value(_liquid_fraction.dot(r, t));
  const auto grad_liquid_fraction = _liquid_fraction.gradient(r, t);
  for (const auto i : make_range(Moose::dim))
    if (_velocity[i])
      liquid_fraction_material_derivative +=
          MetaPhysicL::raw_value((*_velocity[i])(r, t) * grad_liquid_fraction(i));

  return liquid_fraction_material_derivative;
}
