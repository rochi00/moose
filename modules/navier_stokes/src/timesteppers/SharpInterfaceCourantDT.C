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
  params += NonADFunctorInterface::validParams();
  params.addClassDescription(
      "Computes a global timestep from the maximum Courant number of a Rhie-Chow face-flux "
      "provider and, optionally, the accepted liquid-fraction change.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object", "Rhie-Chow user object supplying the current volumetric face flux.");
  params.addRangeCheckedParam<Real>(
      "max_courant", 0.1, "max_courant>0", "Target maximum Courant number.");
  params.addRangeCheckedParam<Real>(
      "initial_dt", 1e-6, "initial_dt>0", "Initial timestep before flux history is available.");
  params.addRangeCheckedParam<Real>(
      "growth_factor", 1.2, "growth_factor>=1", "Maximum timestep growth factor per step.");
  params.addRangeCheckedParam<Real>(
      "max_dt",
      std::numeric_limits<Real>::max(),
      "max_dt>0",
      "Maximum timestep after applying growth and physics-based limits.");
  params.addParam<MooseFunctorName>(
      "liquid_fraction",
      "",
      "Optional liquid-fraction functor used to resolve the diffuse phase-change interval.");
  params.addRangeCheckedParam<Real>(
      "max_liquid_fraction_change",
      0.05,
      "max_liquid_fraction_change>0 & max_liquid_fraction_change<=1",
      "Target maximum accepted cellwise liquid-fraction change.");
  return params;
}

SharpInterfaceCourantDT::SharpInterfaceCourantDT(const InputParameters & parameters)
  : TimeStepper(parameters),
    UserObjectInterface(this),
    NonADFunctorInterface(this),
    _rhie_chow(getUserObject<RhieChowMassFlux>("rhie_chow_user_object")),
    _max_courant(getParam<Real>("max_courant")),
    _initial_dt(getParam<Real>("initial_dt")),
    _growth_factor(getParam<Real>("growth_factor")),
    _max_dt(getParam<Real>("max_dt")),
    _liquid_fraction(isParamValid("liquid_fraction") &&
                             !getParam<MooseFunctorName>("liquid_fraction").empty()
                         ? &getFunctor<Real>("liquid_fraction")
                         : nullptr),
    _max_liquid_fraction_change(getParam<Real>("max_liquid_fraction_change"))
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

  Real phase_change_limited_dt = std::numeric_limits<Real>::max();
  if (_liquid_fraction)
  {
    const auto current_state = Moose::currentState();
    const Moose::StateArg old_state(1, Moose::SolutionIterationType::Time);
    Real maximum_change = 0.0;
    for (const auto * elem_info : _fe_problem.mesh().elemInfoVector())
      if (elem_info->elem()->processor_id() == processor_id())
      {
        const auto elem_arg = makeElemArg(elem_info->elem());
        maximum_change =
            std::max(maximum_change,
                     std::abs((*_liquid_fraction)(elem_arg, current_state) -
                              (*_liquid_fraction)(elem_arg, old_state)));
      }
    _communicator.max(maximum_change);

    if (maximum_change > 0.0)
      phase_change_limited_dt =
          getCurrentDT() * _max_liquid_fraction_change / maximum_change;
  }

  return std::min(
      {_max_dt, courant_limited_dt, phase_change_limited_dt, _growth_factor * getCurrentDT()});
}
