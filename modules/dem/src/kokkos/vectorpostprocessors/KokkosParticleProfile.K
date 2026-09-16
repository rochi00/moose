#include "KokkosParticleProfile.h"
#include "KokkosParticleCloud.h"

#include "MooseMesh.h"

#include "libmesh/mesh_tools.h"

registerMooseObject("DemApp", KokkosParticleProfile);

InputParameters
KokkosParticleProfile::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralVectorPostprocessor::validParams();
  params.addClassDescription(
      "Coarse-grained profiles of a KokkosParticleCloud along a direction: per slab, the number "
      "of particles, the solid fraction, the mean velocity, the granular temperature, and the "
      "stress tensor.");
  params.addRequiredParam<UserObjectName>("cloud", "The KokkosParticleCloud");
  params.addRequiredParam<MooseEnum>("direction", MooseEnum("x y z"), "Direction of the profile");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "num_bins", "num_bins > 0", "Number of slabs, equal in width");
  params.addParam<Real>("lower", "Start of the profile; the mesh bounding box by default");
  params.addParam<Real>("upper", "End of the profile; the mesh bounding box by default");
  return params;
}

KokkosParticleProfile::KokkosParticleProfile(const InputParameters & parameters)
  : Moose::Kokkos::GeneralVectorPostprocessor(parameters),
    _cloud(getUserObject<KokkosParticleCloud>("cloud")),
    _direction(getParam<MooseEnum>("direction")),
    _num_bins(getParam<unsigned int>("num_bins"))
{
  // The slab volume is its width times the bounding box's extent in the other directions of
  // the mesh dimension; particles are spheres of volume 4/3 pi r^3 whatever the dimension
  const auto box = libMesh::MeshTools::create_bounding_box(_fe_problem.mesh().getMesh());
  _lower = isParamValid("lower") ? getParam<Real>("lower") : box.min()(_direction);
  _upper = isParamValid("upper") ? getParam<Real>("upper") : box.max()(_direction);
  if (_upper <= _lower)
    paramError("upper", "Must be above lower");
  _cross_section = 1;
  for (const auto c : make_range(_fe_problem.mesh().dimension()))
    if (c != _direction)
      _cross_section *= box.max()(c) - box.min()(c);
  for (const auto name : {"bin_center",
                          "count",
                          "solid_fraction",
                          "velocity_x",
                          "velocity_y",
                          "velocity_z",
                          "temperature",
                          "stress_xx",
                          "stress_yy",
                          "stress_zz",
                          "stress_xy",
                          "stress_xz",
                          "stress_yz"})
    _columns.push_back(&declareVector(name));
}

void
KokkosParticleProfile::compute()
{
  _cloud.profile(_direction, _lower, _upper, _num_bins, _sums);
}

void
KokkosParticleProfile::finalize()
{
  _communicator.sum(_sums);
  constexpr auto S = KokkosParticleCloud::profile_sums;
  const Real width = (_upper - _lower) / _num_bins;
  const Real volume = width * _cross_section;
  for (auto * column : _columns)
    column->assign(_num_bins, 0);
  for (const auto b : make_range(_num_bins))
  {
    const Real * s = _sums.data() + b * S;
    const Real count = s[0], mass = s[1];
    (*_columns[0])[b] = _lower + (b + 0.5) * width;
    (*_columns[1])[b] = count;
    (*_columns[2])[b] = s[2] / volume;
    Real mean[3] = {0, 0, 0};
    if (mass > 0)
      for (const auto c : make_range(3))
        mean[c] = s[3 + c] / mass;
    for (const auto c : make_range(3))
      (*_columns[3 + c])[b] = mean[c];
    // Fluctuation momentum flux sum m v' v' = sum m v v - mass mean mean, and the temperature
    // m <v'^2> / 3 of the bin's particles
    Real fluct[6];
    constexpr unsigned int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
    for (const auto k : make_range(6))
      fluct[k] = s[6 + k] - mass * mean[ia[k]] * mean[ib[k]];
    (*_columns[6])[b] = count > 0 ? (fluct[0] + fluct[1] + fluct[2]) / (3 * count) : 0;
    for (const auto k : make_range(6))
      (*_columns[7 + k])[b] = -(fluct[k] + s[12 + k]) / volume;
  }
}
