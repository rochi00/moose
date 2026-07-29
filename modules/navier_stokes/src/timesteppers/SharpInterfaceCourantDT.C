#include "SharpInterfaceCourantDT.h"

#include "RhieChowMassFlux.h"

#include <algorithm>
#include <cmath>
#include <limits>

registerMooseObject("NavierStokesApp", SharpInterfaceCourantDT);

InputParameters
SharpInterfaceCourantDT::validParams()
{
  InputParameters params = TimeStepper::validParams();
  params.addClassDescription(
      "Computes a global timestep from the maximum Courant number of a Rhie-Chow face-flux "
      "provider.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object", "Rhie-Chow user object supplying the current volumetric face flux.");
  params.addRangeCheckedParam<Real>(
      "max_courant", 0.1, "max_courant>0", "Target maximum Courant number.");
  params.addRangeCheckedParam<Real>(
      "initial_dt", 1e-6, "initial_dt>0", "Initial timestep before flux history is available.");
  params.addRangeCheckedParam<Real>(
      "growth_factor", 1.2, "growth_factor>=1", "Maximum timestep growth factor per step.");
  return params;
}

SharpInterfaceCourantDT::SharpInterfaceCourantDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    UserObjectInterface(this),
    _rhie_chow(getUserObject<RhieChowMassFlux>("rhie_chow_user_object")),
    _max_courant(getParam<Real>("max_courant")),
    _initial_dt(getParam<Real>("initial_dt")),
    _growth_factor(getParam<Real>("growth_factor"))
{
}

Real
SharpInterfaceCourantDT::computeInitialDT()
{
  return _initial_dt;
}

Real
SharpInterfaceCourantDT::computeDT()
{
  const Real courant_per_unit_dt = _rhie_chow.maxCourant(1.0);
  const Real courant_limited_dt =
      std::isfinite(courant_per_unit_dt) && courant_per_unit_dt > 0.0
          ? _max_courant / courant_per_unit_dt
          : std::numeric_limits<Real>::max();

  return std::min(courant_limited_dt, _growth_factor * getCurrentDT());
}
