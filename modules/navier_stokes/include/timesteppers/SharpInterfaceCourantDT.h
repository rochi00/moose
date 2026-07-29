#pragma once

#include "TimeStepper.h"
#include "UserObjectInterface.h"

class RhieChowMassFlux;

class SharpInterfaceCourantDT : public TimeStepper, public UserObjectInterface
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
};
