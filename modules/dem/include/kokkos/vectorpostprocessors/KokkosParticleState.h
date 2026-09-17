#pragma once

#include "KokkosGeneralVectorPostprocessor.h"

class KokkosParticleCloud;

/**
 * Outputs the position, velocity, and angular velocity of every particle of a
 * KokkosParticleCloud, gathered to the root rank and sorted by global particle ID
 */
class KokkosParticleState : public Moose::Kokkos::GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  KokkosParticleState(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void compute() override;
  virtual void finalize() override;

protected:
  const KokkosParticleCloud & _cloud;

  VectorPostprocessorValue & _gid;
  std::vector<VectorPostprocessorValue *> _x;
  std::vector<VectorPostprocessorValue *> _v;
  std::vector<VectorPostprocessorValue *> _omega;
  /// Force and torque, only with output_forces
  std::vector<VectorPostprocessorValue *> _f;
  std::vector<VectorPostprocessorValue *> _tau;
  /// Radius, only with output_radius
  VectorPostprocessorValue * _r = nullptr;
};
