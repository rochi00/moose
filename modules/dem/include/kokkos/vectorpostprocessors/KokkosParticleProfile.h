#pragma once

#include "KokkosGeneralVectorPostprocessor.h"

class KokkosParticleCloud;

/**
 * Coarse-grained profiles of a KokkosParticleCloud along one direction: per slab, the number of
 * particles, the solid fraction, the mean velocity, the granular temperature, and the stress
 * tensor (kinetic plus contact), all reduced across ranks
 */
class KokkosParticleProfile : public Moose::Kokkos::GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  KokkosParticleProfile(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void compute() override;
  virtual void finalize() override;

protected:
  const KokkosParticleCloud & _cloud;
  const unsigned int _direction;
  const std::size_t _num_bins;
  /// Range of the profile and the cross-section of a slab (the bounding box's other extents)
  Real _lower, _upper, _cross_section;
  /// Local sums from the cloud, reduced in finalize()
  std::vector<Real> _sums;
  /// Output columns: bin center, count, solid fraction, velocity (3), temperature, stress (6)
  std::vector<VectorPostprocessorValue *> _columns;
};
