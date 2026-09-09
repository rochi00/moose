//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearWCNSFVDispersePhaseDragFunctorMaterial.h"
#include "NavierStokesMethods.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", LinearWCNSFVDispersePhaseDragFunctorMaterial);

InputParameters
LinearWCNSFVDispersePhaseDragFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Computes the linear drag function between the dispersed and the continuous phase, from the "
      "particle Reynolds number formed with the slip velocity, for the linear finite volume "
      "discretization.");
  params.addParam<MooseFunctorName>("drag_coef_name",
                                    "Darcy_coefficient",
                                    "Name of the scalar friction coefficient defined. The vector "
                                    "coefficient is suffixed with _vec");
  params.addRequiredParam<MooseFunctorName>("u", "The slip velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The slip velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The slip velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(
      NS::density,
      "Continuous phase density, which is the density the particle Reynolds number is formed "
      "from.");
  params.addRequiredParam<MooseFunctorName>(
      NS::mu,
      "Continuous phase dynamic viscosity, which is the viscosity the particle Reynolds number is "
      "formed from.");
  params.addParam<MooseFunctorName>(
      "particle_diameter", 1.0, "Diameter of particles in the dispersed phase.");
  return params;
}

LinearWCNSFVDispersePhaseDragFunctorMaterial::LinearWCNSFVDispersePhaseDragFunctorMaterial(
    const InputParameters & parameters)
  : FunctorMaterial(parameters),
    NonADFunctorInterface(this),
    _dim(_subproblem.mesh().dimension()),
    _u_var(NonADFunctorInterface::getFunctor<Real>("u")),
    _v_var(parameters.isParamValid("v") ? &NonADFunctorInterface::getFunctor<Real>("v") : nullptr),
    _w_var(parameters.isParamValid("w") ? &NonADFunctorInterface::getFunctor<Real>("w") : nullptr),
    _rho(NonADFunctorInterface::getFunctor<Real>(NS::density)),
    _mu(NonADFunctorInterface::getFunctor<Real>(NS::mu)),
    _particle_diameter(NonADFunctorInterface::getFunctor<Real>("particle_diameter"))
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied.");
  if (_dim >= 3 && !_w_var)
    paramError("w", "In three dimensions, the w velocity must be supplied.");

  const auto drag = [this](const auto & r, const auto & t) -> Real
  {
    // The particle Reynolds number is formed from the slip velocity and the continuous phase
    // properties, Re_p = rho_c d_d |u_slip| / mu_c. See Manninen, Taivassalo and Kallio, VTT
    // Publications 288 (1996), equation (39).
    //
    // Reading the slip velocity here does not close a loop back onto this object:
    // LinearWCNSFV2PSlipVelocityFunctorMaterial solves the drag correlation together with its
    // own force balance rather than consuming a drag functor.
    RealVectorValue slip_velocity(_u_var(r, t));
    if (_dim > 1)
      slip_velocity(1) = (*_v_var)(r, t);
    if (_dim > 2)
      slip_velocity(2) = (*_w_var)(r, t);
    const auto slip_speed = NS::computeSpeed<Real>(slip_velocity);

    return NS::dragFunction(
        NS::particleReynoldsNumber(_rho(r, t), _particle_diameter(r, t), slip_speed, _mu(r, t)));
  };

  const auto & drag_func =
      addFunctorProperty<Real>(getParam<MooseFunctorName>("drag_coef_name"), drag);

  // The vector form is what the momentum friction kernels consume
  const auto drag_vec = [&drag_func](const auto & r, const auto & t) -> RealVectorValue
  {
    const auto value = drag_func(r, t);
    return RealVectorValue(value, value, value);
  };
  addFunctorProperty<RealVectorValue>(getParam<MooseFunctorName>("drag_coef_name") + "_vec",
                                      drag_vec);
}
