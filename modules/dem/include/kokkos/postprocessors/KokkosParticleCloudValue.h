#pragma once

#include "KokkosGeneralPostprocessor.h"

class KokkosParticleCloud;

/**
 * Reports one of the global scalars of a KokkosParticleCloud
 */
class KokkosParticleCloudValue : public Moose::Kokkos::GeneralPostprocessor
{
public:
  static InputParameters validParams();

  KokkosParticleCloudValue(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void compute() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  const KokkosParticleCloud & _cloud;
  const MooseEnum _value;
};
