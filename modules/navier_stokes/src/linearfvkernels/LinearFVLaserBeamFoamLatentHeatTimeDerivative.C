//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVLaserBeamFoamLatentHeatTimeDerivative.h"

registerMooseObject("NavierStokesApp", LinearFVLaserBeamFoamLatentHeatTimeDerivative);
registerMooseObjectAliased("NavierStokesApp",
                           LinearFVLaserBeamFoamLatentHeatTimeDerivative,
                           "LinearFVEnthalpyLiquidFractionLatentHeatTimeDerivative");

InputParameters
LinearFVLaserBeamFoamLatentHeatTimeDerivative::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Adds the LaserBeamFoam-style latent heat RHS contribution "
      "-L * d(rho * liquid_fraction) / dt.");
  params.addRequiredParam<MooseFunctorName>("L", "Latent heat.");
  params.addRequiredParam<MooseFunctorName>("density", "Density multiplier.");
  params.addParam<MooseFunctorName>(
      "density_dot",
      "Optional density time derivative functor. When supplied, the kernel evaluates "
      "d(rho * liquid_fraction) / dt = rho * liquid_fraction_dot + density_dot * "
      "liquid_fraction.");
  params.addRequiredParam<MooseFunctorName>("liquid_fraction", "Liquid fraction.");
  return params;
}

LinearFVLaserBeamFoamLatentHeatTimeDerivative::
    LinearFVLaserBeamFoamLatentHeatTimeDerivative(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _L(getFunctor<Real>("L")),
    _density(getFunctor<Real>("density")),
    _density_dot(isParamValid("density_dot") ? &getFunctor<Real>("density_dot") : nullptr),
    _liquid_fraction(getFunctor<Real>("liquid_fraction"))
{
}

Real
LinearFVLaserBeamFoamLatentHeatTimeDerivative::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVLaserBeamFoamLatentHeatTimeDerivative::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto state = determineState();
  const Real density = _density(elem_arg, state);
  const Real liquid_fraction = _liquid_fraction(elem_arg, state);
  const Real density_liquid_fraction_dot =
      density * _liquid_fraction.dot(elem_arg, state) +
      (_density_dot ? (*_density_dot)(elem_arg, state) * liquid_fraction : 0.0);

  return -_L(elem_arg, state) * density_liquid_fraction_dot * _current_elem_volume;
}
