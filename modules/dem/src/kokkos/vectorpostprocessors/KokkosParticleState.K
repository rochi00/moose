#include "KokkosParticleState.h"
#include "KokkosParticleCloud.h"

registerMooseObject("DemApp", KokkosParticleState);

InputParameters
KokkosParticleState::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralVectorPostprocessor::validParams();
  params.addClassDescription("Outputs the position, velocity, and angular velocity of every "
                             "particle of a KokkosParticleCloud, sorted by global particle ID.");
  params.addRequiredParam<UserObjectName>("cloud", "The KokkosParticleCloud to output");
  params.addParam<bool>(
      "output_forces", false, "Also output the force and body-frame torque on every particle");
  return params;
}

KokkosParticleState::KokkosParticleState(const InputParameters & parameters)
  : Moose::Kokkos::GeneralVectorPostprocessor(parameters),
    _cloud(getUserObject<KokkosParticleCloud>("cloud")),
    _gid(declareVector("gid")),
    _x({&declareVector("x"), &declareVector("y"), &declareVector("z")}),
    _v({&declareVector("vx"), &declareVector("vy"), &declareVector("vz")}),
    _omega({&declareVector("wx"), &declareVector("wy"), &declareVector("wz")})
{
  if (getParam<bool>("output_forces"))
  {
    _f = {&declareVector("fx"), &declareVector("fy"), &declareVector("fz")};
    _tau = {&declareVector("taux"), &declareVector("tauy"), &declareVector("tauz")};
  }
}

void
KokkosParticleState::compute()
{
  const auto & cloud = _cloud.cloud();
  const auto n = cloud.size();
  auto gid = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cloud.gid);
  auto x = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cloud.x);
  auto v = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cloud.v);
  auto omega = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cloud.omega);
  auto f = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cloud.f_out);
  auto tau = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, cloud.tau_out);

  _gid.resize(n);
  for (const auto c : make_range(3))
  {
    _x[c]->resize(n);
    _v[c]->resize(n);
    _omega[c]->resize(n);
    if (_f.size())
    {
      _f[c]->resize(n);
      _tau[c]->resize(n);
    }
  }
  for (const auto i : make_range(n))
  {
    _gid[i] = gid(i);
    for (const auto c : make_range(3))
    {
      (*_x[c])[i] = x(i, c);
      (*_v[c])[i] = v(i, c);
      (*_omega[c])[i] = omega(i, c);
      if (_f.size())
      {
        (*_f[c])[i] = f(i, c);
        (*_tau[c])[i] = tau(i, c);
      }
    }
  }
}

void
KokkosParticleState::finalize()
{
  _communicator.gather(0, _gid);
  for (const auto c : make_range(3))
  {
    _communicator.gather(0, *_x[c]);
    _communicator.gather(0, *_v[c]);
    _communicator.gather(0, *_omega[c]);
    if (_f.size())
    {
      _communicator.gather(0, *_f[c]);
      _communicator.gather(0, *_tau[c]);
    }
  }

  // Sort by global ID so the output does not depend on the partitioning
  std::vector<std::size_t> order(_gid.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [this](auto a, auto b) { return _gid[a] < _gid[b]; });
  auto permute = [&order](VectorPostprocessorValue & values)
  {
    VectorPostprocessorValue sorted(values.size());
    for (const auto i : index_range(order))
      sorted[i] = values[order[i]];
    values.swap(sorted);
  };
  permute(_gid);
  for (const auto c : make_range(3))
  {
    permute(*_x[c]);
    permute(*_v[c]);
    permute(*_omega[c]);
    if (_f.size())
    {
      permute(*_f[c]);
      permute(*_tau[c]);
    }
  }
}
