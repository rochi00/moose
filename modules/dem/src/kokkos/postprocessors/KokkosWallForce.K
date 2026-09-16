#include "KokkosWallForce.h"
#include "KokkosParticleCloud.h"

registerMooseObject("DemApp", KokkosWallForce);

InputParameters
KokkosWallForce::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralPostprocessor::validParams();
  params.addClassDescription("Reports the total contact force the particles of a "
                             "KokkosParticleCloud exert on one of its walls.");
  params.addRequiredParam<UserObjectName>("cloud", "The KokkosParticleCloud");
  params.addParam<unsigned int>("wall", "Index of an analytic wall (in wall_points order)");
  params.addParam<BoundaryName>("boundary", "A sideset wall (one of wall_boundaries)");
  params.addParam<MooseEnum>("component",
                             MooseEnum("x y z magnitude", "magnitude"),
                             "Component of the force, or its magnitude");
  return params;
}

KokkosWallForce::KokkosWallForce(const InputParameters & parameters)
  : Moose::Kokkos::GeneralPostprocessor(parameters),
    _cloud(getUserObject<KokkosParticleCloud>("cloud")),
    _component(getParam<MooseEnum>("component"))
{
  if (isParamValid("wall") == isParamValid("boundary"))
    paramError("wall", "Give either an analytic wall index or a sideset boundary");
  // The cloud's walls exist once it is constructed; both are constructed before any executes
  if (isParamValid("wall"))
    _wall = getParam<unsigned int>("wall");
  else
    _wall = _cloud.wallIndex(getParam<BoundaryName>("boundary"));
  if (_wall >= _cloud.numWalls())
    paramError("wall", "The cloud has ", _cloud.numWalls(), " walls");
}

PostprocessorValue
KokkosWallForce::getValue() const
{
  const auto force = _cloud.wallForce(_wall);
  if (_component == "magnitude")
    return force.norm();
  return force(_component);
}
