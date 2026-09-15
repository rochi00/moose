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
                                     MooseEnum("num_particles num_exited num_unresolved num_migrated "
                                               "max_hops num_neighbor_pairs num_neighbor_list_builds "
                                               "num_ghosts num_ghost_forwards num_contacts "
                                               "coordination_number load_imbalance kinetic_energy "
                                               "rotational_kinetic_energy potential_energy "
                                               "total_energy linear_momentum angular_momentum "
                                               "virial_pressure"),
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
  if (_value == "num_exited")
    return _cloud.numExited();
  if (_value == "num_unresolved")
    return _cloud.numUnresolved();
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
  if (_value == "num_ghost_forwards")
    return _cloud.numGhostForwards();
  if (_value == "num_contacts")
    return _cloud.numContacts();
  if (_value == "coordination_number")
    return _cloud.coordinationNumber();
  if (_value == "load_imbalance")
    return _cloud.loadImbalance();
  if (_value == "kinetic_energy")
    return _cloud.kineticEnergy();
  if (_value == "rotational_kinetic_energy")
    return _cloud.rotationalKineticEnergy();
  if (_value == "potential_energy")
    return _cloud.potentialEnergy();
  if (_value == "total_energy")
    return _cloud.totalEnergy();
  if (_value == "linear_momentum")
    return _cloud.linearMomentum();
  if (_value == "angular_momentum")
    return _cloud.angularMomentum();
  return _cloud.virialPressure();
}
