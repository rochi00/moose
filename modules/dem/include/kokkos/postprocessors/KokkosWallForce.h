#pragma once

#include "KokkosGeneralPostprocessor.h"

class KokkosParticleCloud;

/**
 * A component or the magnitude of the total contact force the particles of a KokkosParticleCloud
 * exert on one of its walls: an analytic wall by index, or a sideset wall by boundary name
 */
class KokkosWallForce : public Moose::Kokkos::GeneralPostprocessor
{
public:
  static InputParameters validParams();

  KokkosWallForce(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void compute() override {}
  virtual void finalize() override {}
  virtual PostprocessorValue getValue() const override;

protected:
  const KokkosParticleCloud & _cloud;
  /// Index of the wall in the cloud's numbering, resolved at setup
  std::size_t _wall = 0;
  const MooseEnum _component;
};
