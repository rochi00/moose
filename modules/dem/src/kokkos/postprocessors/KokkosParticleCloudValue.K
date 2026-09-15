#include "KokkosParticleCloudValue.h"
#include "KokkosParticleCloud.h"

registerMooseObject("DemApp", KokkosParticleCloudValue);

InputParameters
KokkosParticleCloudValue::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralPostprocessor::validParams();
  params.addClassDescription("Reports one of the global scalars of a KokkosParticleCloud.");
  params.addRequiredParam<UserObjectName>("cloud", "The KokkosParticleCloud to report on");
  params.addRequiredParam<MooseEnum>("value",
                                     MooseEnum("num_particles num_lost num_migrated max_hops "
                                               "num_neighbor_pairs num_neighbor_list_builds num_ghosts "
                                               "kinetic_energy rotational_kinetic_energy "
                                               "angular_momentum"),
                                     "Which value to report");
  return params;
}

KokkosParticleCloudValue::KokkosParticleCloudValue(const InputParameters & parameters)
  : Moose::Kokkos::GeneralPostprocessor(parameters),
    _cloud(getUserObject<KokkosParticleCloud>("cloud")),
    _value(getParam<MooseEnum>("value"))
{
}

PostprocessorValue
KokkosParticleCloudValue::getValue() const
{
  if (_value == "num_particles")
    return _cloud.numParticles();
  if (_value == "num_lost")
    return _cloud.numLost();
  if (_value == "num_migrated")
    return _cloud.numMigrated();
  if (_value == "max_hops")
    return _cloud.maxHops();
  if (_value == "num_neighbor_pairs")
    return _cloud.numNeighborPairs();
  if (_value == "num_neighbor_list_builds")
    return _cloud.numNeighborListBuilds();
  if (_value == "num_ghosts")
    return _cloud.numGhosts();
  if (_value == "kinetic_energy")
    return _cloud.kineticEnergy();
  if (_value == "rotational_kinetic_energy")
    return _cloud.rotationalKineticEnergy();
  return _cloud.angularMomentum();
}
