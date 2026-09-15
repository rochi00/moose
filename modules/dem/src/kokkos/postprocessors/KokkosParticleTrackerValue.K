#include "KokkosParticleTrackerValue.h"
#include "KokkosParticleTracker.h"

registerMooseObject("DemApp", KokkosParticleTrackerValue);

InputParameters
KokkosParticleTrackerValue::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralPostprocessor::validParams();
  params.addClassDescription("Reports one of the counters of a KokkosParticleTracker.");
  params.addRequiredParam<UserObjectName>("tracker", "The KokkosParticleTracker to report on");
  params.addRequiredParam<MooseEnum>(
      "value",
      MooseEnum("num_particles num_lost num_migrated max_hops"),
      "Which counter to report");
  return params;
}

KokkosParticleTrackerValue::KokkosParticleTrackerValue(const InputParameters & parameters)
  : Moose::Kokkos::GeneralPostprocessor(parameters),
    _tracker(getUserObject<KokkosParticleTracker>("tracker")),
    _value(getParam<MooseEnum>("value"))
{
}

PostprocessorValue
KokkosParticleTrackerValue::getValue() const
{
  if (_value == "num_particles")
    return _tracker.numParticles();
  if (_value == "num_lost")
    return _tracker.numLost();
  if (_value == "num_migrated")
    return _tracker.numMigrated();
  return _tracker.maxHops();
}
