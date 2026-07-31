#pragma once

#include "NonADFunctorInterface.h"
#include "TimeStepper.h"
#include "UserObjectInterface.h"

class RhieChowMassFlux;

class SharpInterfaceCourantDT : public TimeStepper,
                                public UserObjectInterface,
                                public NonADFunctorInterface
{
public:
  static InputParameters validParams();

  SharpInterfaceCourantDT(const InputParameters & parameters);

protected:
  virtual Real computeInitialDT() override;
  virtual Real computeDT() override;

  const RhieChowMassFlux & _rhie_chow;
  const Real _max_courant;
  const Real _initial_dt;
  const Real _growth_factor;
  const Real _max_dt;
  const Moose::Functor<Real> * const _liquid_fraction;
  const Real _max_liquid_fraction_change;
};
