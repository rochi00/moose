#pragma once

#include "KokkosGeneralPostprocessor.h"

class KokkosParticleTracker;

/**
 * Reports one of the counters of a KokkosParticleTracker
 */
class KokkosParticleTrackerValue : public Moose::Kokkos::GeneralPostprocessor
{
public:
  static InputParameters validParams();

  KokkosParticleTrackerValue(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void compute() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  const KokkosParticleTracker & _tracker;
  const MooseEnum _value;
};
