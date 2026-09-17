#include "KokkosParticleCloud.h"

#include "MooseMesh.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_tools.h"
#include "timpi/parallel_sync.h"

#include <algorithm>

registerMooseObject("DemApp", KokkosParticleCloud);

/// Reals forwarded per ghost every substep: position, velocity, and spin
static constexpr std::size_t ghost_stride = 9;
/// Reals a contact history migrates as: the two key IDs and the tangential and rolling springs
static constexpr std::size_t pair_state_record_size = 8;

InputParameters
KokkosParticleCloud::validParams()
{
  InputParameters params = Moose::Kokkos::GeneralUserObject::validParams();
  params.addClassDescription(
      "Integrates a cloud of spherical particles on device with Velocity Verlet over a number of "
      "substeps per time step, tracking the element containing each particle by walking across "
      "element faces and migrating particles between ranks as they cross partitions.");
  params.addRequiredParam<std::vector<Point>>("initial_positions",
                                              "Initial positions of all particles");
  params.addParam<RealVectorValue>(
      "initial_velocity", RealVectorValue(0, 0, 0), "Initial velocity of every particle");
  params.addParam<std::vector<Point>>(
      "initial_velocities",
      "Per-particle initial velocities overriding initial_velocity, one per initial position");
  params.addParam<std::vector<unsigned int>>(
      "initial_frozen",
      {},
      "Indices (into initial_positions) of particles that are frozen in place: never moved, "
      "immovable obstacles to the others (LAMMPS fix freeze)");
  params.addParam<RealVectorValue>("initial_angular_velocity",
                                   RealVectorValue(0, 0, 0),
                                   "Initial body-frame angular velocity of every particle");
  params.addRequiredRangeCheckedParam<Real>("radius", "radius > 0", "Particle radius");
  params.addParam<std::vector<Real>>(
      "initial_radii", "Per-particle radii overriding radius, one per initial position");
  params.addRequiredRangeCheckedParam<Real>("density", "density > 0", "Particle density");
  params.addParam<bool>("disk_mass",
                        false,
                        "Give the particles the mass of disks, density times pi r^2, as LAMMPS "
                        "does in two dimensions (the moment of inertia stays 0.4 m r^2, as "
                        "there); otherwise the mass of spheres");
  params.addParam<RealVectorValue>(
      "gravity", RealVectorValue(0, 0, 0), "Gravitational acceleration");
  params.addParam<MooseEnum>("contact_model",
                             MooseEnum("linear_spring_dashpot hertz", "linear_spring_dashpot"),
                             "Normal contact model between spheres and against walls");
  params.addRangeCheckedParam<Real>("normal_stiffness",
                                    0,
                                    "normal_stiffness >= 0",
                                    "Spring stiffness of the linear spring-dashpot normal contact; "
                                    "zero disables contact");
  params.addRangeCheckedParam<Real>("normal_damping",
                                    0,
                                    "normal_damping >= 0",
                                    "Dashpot coefficient of the linear spring-dashpot normal contact");
  params.addRangeCheckedParam<Real>("tangential_stiffness",
                                    0,
                                    "tangential_stiffness >= 0",
                                    "Tangential spring stiffness of the linear spring-dashpot "
                                    "contact, zero for frictionless contact; with the Hertz "
                                    "model, when given, the coefficient 8 G* of the Mindlin "
                                    "stiffness 8 G* sqrt(r_eff overlap) in place of the one "
                                    "from the elastic constants (LAMMPS's k_t)");
  params.addRangeCheckedParam<Real>("tangential_damping",
                                    0,
                                    "tangential_damping >= 0",
                                    "Tangential dashpot coefficient of the linear spring-dashpot "
                                    "contact");
  params.addParam<bool>("rescale_histories",
                        true,
                        "Keep the length of the tangential and rolling springs when they are "
                        "rotated into the current tangent plane of a turning contact (LAMMPS "
                        "pair granular); otherwise only project them (LAMMPS gran/* styles, "
                        "LIGGGHTS)");
  params.addParam<bool>("deformed_torque_arm",
                        false,
                        "Apply the tangential force at the contact plane, r_i - overlap / 2 "
                        "from the center, when computing the torque (LAMMPS pair granular); "
                        "otherwise at the undeformed radius");
  params.addParam<bool>("mass_scaled_damping",
                        false,
                        "Multiply the linear model's normal_damping and tangential_damping by "
                        "the effective mass of the contact, as LAMMPS's gran/* styles do with "
                        "their gamma (a rate) rather than a dashpot coefficient");
  params.addParam<bool>("limit_damping",
                        false,
                        "Never let the normal damping make the normal force attractive: the "
                        "force is clamped at zero while the spring's repulsion is smaller than "
                        "the dashpot's pull, as in LAMMPS's granular styles with limit_damping");
  params.addRangeCheckedParam<Real>(
      "friction",
      0,
      "friction >= 0",
      "Coefficient of friction limiting the tangential force to this times the normal force; "
      "zero for frictionless contact");
  params.addParam<bool>("partial_slip",
                        false,
                        "Follow the Mindlin-Deresiewicz partial-slip loading curve for the "
                        "tangential spring, which reaches the friction limit at 3/2 of the "
                        "no-slip slip displacement, instead of the linear no-slip spring");
  params.addRangeCheckedParam<Real>("rolling_friction",
                                    0,
                                    "rolling_friction >= 0",
                                    "Coefficient of rolling friction limiting the rolling "
                                    "resistance torque to this times the normal force times the "
                                    "effective radius; zero for no rolling resistance");
  params.addRangeCheckedParam<Real>("rolling_stiffness",
                                    0,
                                    "rolling_stiffness >= 0",
                                    "Rolling spring stiffness, on the rolling displacement");
  params.addRangeCheckedParam<Real>(
      "rolling_damping", 0, "rolling_damping >= 0", "Rolling dashpot coefficient");
  params.addRangeCheckedParam<Real>("wall_friction",
                                    "wall_friction >= 0",
                                    "Coefficient of friction of the particle-wall contacts; "
                                    "defaults to friction");
  params.addRangeCheckedParam<Real>("wall_rolling_friction",
                                    "wall_rolling_friction >= 0",
                                    "Coefficient of rolling friction of the particle-wall "
                                    "contacts; defaults to rolling_friction");
  params.addRangeCheckedParam<Real>(
      "youngs_modulus", "youngs_modulus > 0", "Young's modulus of the Hertz normal contact");
  params.addRangeCheckedParam<Real>("poissons_ratio",
                                    "poissons_ratio >= 0 & poissons_ratio < 0.5",
                                    "Poisson's ratio of the Hertz normal contact");
  params.addRangeCheckedParam<Real>("restitution",
                                    1,
                                    "restitution > 0 & restitution <= 1",
                                    "Coefficient of restitution the Tsuji damping of the Hertz "
                                    "normal contact is set to reproduce; one for no damping");
  params.addRangeCheckedParam<Real>(
      "dissipation_time",
      0,
      "dissipation_time >= 0",
      "Dissipative constant A of the Kuwabara-Kono viscoelastic damping of the Hertz normal "
      "contact, used instead of the Tsuji damping when positive: the damping force is A times "
      "the rate of the elastic force, and the restitution coefficient falls with the impact "
      "velocity as 1 - e ~ v^(1/5)");
  params.addParam<std::vector<RealVectorValue>>(
      "wall_velocities", {}, "Translation velocity of each fixed planar wall, one per wall");
  params.addParam<Point>("wall_rotation_center",
                         Point(0, 0, 0),
                         "Point on the axis about which the analytic walls rotate rigidly");
  params.addParam<RealVectorValue>("wall_angular_velocity",
                                   RealVectorValue(0, 0, 0),
                                   "Angular velocity (axis times rate, in radians per unit "
                                   "time) of the rigid rotation of all the analytic walls about "
                                   "wall_rotation_center; the walls' points and normals turn "
                                   "with it and their surface velocity includes it");
  params.addParam<std::vector<unsigned int>>(
      "servo_walls", {}, "Planar walls (by index) driven to a target normal force by a servo");
  params.addParam<std::vector<Real>>(
      "servo_forces",
      {},
      "Target compressive force on each servo wall (the particles' push against its normal)");
  params.addRangeCheckedParam<Real>(
      "servo_gain",
      0,
      "servo_gain >= 0",
      "Servo velocity per unit of force error (m/s per N): the wall moves into the domain when "
      "the force is below the target and out when above");
  params.addRangeCheckedParam<Real>(
      "servo_max_velocity", 0, "servo_max_velocity >= 0", "Largest servo velocity");
  params.addParam<std::vector<BoundaryName>>(
      "wall_boundaries",
      {},
      "Sidesets whose faces are walls, triangulated, with the mesh interior as the inside");
  params.addParam<std::vector<VariableName>>(
      "wall_displacements",
      {},
      "Nodal displacement variables (one per mesh dimension) moving the sideset walls: every "
      "MOOSE step the faces are carried to the displaced positions at the velocity that gets "
      "them there over the step");
  params.addRangeCheckedParam<Real>(
      "wall_reach",
      0,
      "wall_reach >= 0",
      "Extra distance within which a sideset face is a candidate of an element, at least the "
      "largest wall displacement expected");
  params.addRangeCheckedParam<Real>(
      "wall_curvature",
      1,
      "wall_curvature >= 0 & wall_curvature < 180",
      "Largest angle in degrees between the normals of two sideset faces sharing an edge for "
      "them to be one surface, across which a contact history is carried as the contact point "
      "slides from one to the other");
  params.addParam<std::vector<Point>>(
      "wall_points", {}, "A point on each fixed planar wall, one per entry of wall_normals");
  params.addParam<std::vector<RealVectorValue>>(
      "wall_normals",
      {},
      "Normal of each fixed planar wall pointing into the domain, one per entry of wall_points; "
      "particles overlapping a wall get the normal contact force against it");
  params.addParam<MultiMooseEnum>(
      "periodic",
      MultiMooseEnum("x y z"),
      "Directions in which the domain is periodic over the extent of the mesh bounding box");
  params.addRangeCheckedParam<unsigned int>(
      "substeps", 1, "substeps > 0", "Number of integration substeps per time step");
  // A particle should cross only a few cells per substep; a much larger walk indicates a mesh the
  // face-plane test cannot handle, so bound it rather than loop for the mesh size
  params.addParam<unsigned int>(
      "max_hops", 32, "Maximum number of face-neighbor hops per particle per substep");
  params.addParam<MooseEnum>("broad_phase",
                             MooseEnum("uniform_grid arborx_bvh", "uniform_grid"),
                             "Broad phase of the neighbor list: a uniform grid of cell size "
                             "2 r_max + skin, or an ArborX bounding volume hierarchy whose cost "
                             "does not depend on the size distribution (needs the ArborX "
                             "submodule)");
  params.addRangeCheckedParam<Real>("skin",
                                    "skin >= 0",
                                    "Verlet skin distance of the neighbor list; the list is rebuilt "
                                    "once some particle has moved skin / 2. Defaults to a tenth of "
                                    "radius");
  params.addParam<MooseEnum>("timestep_check",
                             MooseEnum("none warn error", "warn"),
                             "Check that the substep is below timestep_fraction of the contact "
                             "time scales (the contact duration of the linear model, or the "
                             "Rayleigh and Hertz times of the Hertz model, as LIGGGHTS's "
                             "check/timestep/gran) at setup and every step, and warn or error "
                             "when it is not");
  params.addRangeCheckedParam<Real>(
      "timestep_fraction",
      0.2,
      "timestep_fraction > 0",
      "Largest allowed ratio of the substep to the contact time scales");
  params.addParam<std::vector<Point>>(
      "insertion_box",
      {},
      "Two corners of a box in which particles are inserted at insertion_rate, at random "
      "positions overlapping no particle, with the given radius and density");
  params.addRangeCheckedParam<Real>(
      "insertion_rate", 0, "insertion_rate >= 0", "Particles inserted per unit time");
  params.addParam<std::vector<Real>>(
      "insertion_radius_range",
      {},
      "Lower and upper radius of the inserted particles, drawn uniformly (LAMMPS fix pour "
      "diam range); the upper one may not exceed radius. Every inserted particle has radius "
      "when not given");
  params.addParam<RealVectorValue>(
      "insertion_velocity", RealVectorValue(0, 0, 0), "Velocity of the inserted particles");
  params.addParam<Real>("insertion_start_time", 0, "Time insertion starts");
  params.addParam<Real>("insertion_end_time",
                        std::numeric_limits<Real>::max(),
                        "Time insertion ends");
  params.addParam<unsigned int>("insertion_seed", 0, "Seed of the insertion positions");
  params.addParam<std::vector<Point>>(
      "outflow_box", {}, "Two corners of a box in which particles are removed");
  params.addParam<bool>("allow_unresolved",
                        false,
                        "Allow particles to remain unresolved (not located in any element) at "
                        "the end of a step instead of erroring");
  params.addParam<bool>("gpu_aware_mpi",
                        false,
                        "Pass device buffers directly to MPI for the ghost update on a GPU build; "
                        "requires a GPU-aware MPI library");
  params.addParam<bool>("verify",
                        false,
                        "Check the device element assignment against libMesh's PointLocator "
                        "every step and error on any mismatch");
  params.addParam<bool>("verify_neighbor_list",
                        false,
                        "Check the neighbor list against a brute-force pair search after every "
                        "build and error on any mismatch");
  params.addParam<bool>("verify_pair_states",
                        false,
                        "Check every step that the contact histories with a stretched tangential "
                        "spring are exactly those of the overlapping pairs, and error otherwise");
  return params;
}

KokkosParticleCloud::KokkosParticleCloud(const InputParameters & parameters)
  : Moose::Kokkos::GeneralUserObject(parameters),
    _initial_positions(getParam<std::vector<Point>>("initial_positions")),
    _initial_velocity(getParam<RealVectorValue>("initial_velocity")),
    _initial_velocities(isParamValid("initial_velocities")
                            ? getParam<std::vector<Point>>("initial_velocities")
                            : std::vector<Point>()),
    _initial_angular_velocity(getParam<RealVectorValue>("initial_angular_velocity")),
    _radius(getParam<Real>("radius")),
    _initial_radii(isParamValid("initial_radii") ? getParam<std::vector<Real>>("initial_radii")
                                                 : std::vector<Real>()),
    _initial_frozen(getParam<std::vector<unsigned int>>("initial_frozen").begin(),
                    getParam<std::vector<unsigned int>>("initial_frozen").end()),
    _density(getParam<Real>("density")),
    _disk_mass(getParam<bool>("disk_mass")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _contact_model(getParam<MooseEnum>("contact_model").getEnum<ContactModel>()),
    _linear({getParam<Real>("normal_stiffness"),
             getParam<Real>("normal_damping"),
             getParam<Real>("tangential_stiffness"),
             getParam<Real>("tangential_damping"),
             getParam<bool>("limit_damping"),
             getParam<bool>("mass_scaled_damping")}),
    _hertz(makeHertz()),
    _friction({getParam<Real>("friction"),
               getParam<bool>("partial_slip"),
               getParam<Real>("rolling_friction"),
               getParam<Real>("rolling_stiffness"),
               getParam<Real>("rolling_damping"),
               getParam<bool>("rescale_histories")}),
    _wall_friction({isParamValid("wall_friction") ? getParam<Real>("wall_friction")
                                                  : getParam<Real>("friction"),
                    getParam<bool>("partial_slip"),
                    isParamValid("wall_rolling_friction") ? getParam<Real>("wall_rolling_friction")
                                                          : getParam<Real>("rolling_friction"),
                    getParam<Real>("rolling_stiffness"),
                    getParam<Real>("rolling_damping"),
                    getParam<bool>("rescale_histories")}),
    _deformed_torque_arm(getParam<bool>("deformed_torque_arm")),
    _walls(getParam<std::vector<Point>>("wall_points"),
           getParam<std::vector<RealVectorValue>>("wall_normals")),
    _wall_velocities(getParam<std::vector<RealVectorValue>>("wall_velocities")),
    _servo_walls(getParam<std::vector<unsigned int>>("servo_walls")),
    _servo_forces(getParam<std::vector<Real>>("servo_forces")),
    _servo_gain(getParam<Real>("servo_gain")),
    _servo_max_velocity(getParam<Real>("servo_max_velocity")),
    _wall_boundaries(getParam<std::vector<BoundaryName>>("wall_boundaries")),
    _wall_reach(getParam<Real>("wall_reach")),
    _wall_curvature(getParam<Real>("wall_curvature") * libMesh::pi / 180),
    _periodic(getParam<MultiMooseEnum>("periodic")),
    _substeps(getParam<unsigned int>("substeps")),
    _max_hops(getParam<unsigned int>("max_hops")),
    _skin(isParamValid("skin") ? getParam<Real>("skin") : 0.1 * _radius),
    _broad_phase(getParam<MooseEnum>("broad_phase").getEnum<DEM::BroadPhase>()),
    _verify(getParam<bool>("verify")),
    _verify_neighbor_list(getParam<bool>("verify_neighbor_list")),
    _verify_pair_states(getParam<bool>("verify_pair_states")),
    _timestep_check(getParam<MooseEnum>("timestep_check")),
    _timestep_fraction(getParam<Real>("timestep_fraction")),
    _allow_unresolved(getParam<bool>("allow_unresolved")),
    _insertion_box(getParam<std::vector<Point>>("insertion_box")),
    _insertion_rate(getParam<Real>("insertion_rate")),
    _insertion_radius_range(getParam<std::vector<Real>>("insertion_radius_range")),
    _insertion_velocity(getParam<RealVectorValue>("insertion_velocity")),
    _insertion_start_time(getParam<Real>("insertion_start_time")),
    _insertion_end_time(getParam<Real>("insertion_end_time")),
    _insertion_generator(getParam<unsigned int>("insertion_seed")),
    _outflow_box(getParam<std::vector<Point>>("outflow_box")),
    _checkpoint(declareRestartableData<DEM::CloudCheckpoint>("checkpoint")),
    _ghost_tag(_communicator.get_unique_tag()),
    _gpu_aware_mpi(getParam<bool>("gpu_aware_mpi"))
{
  if (_initial_radii.size() && _initial_radii.size() != _initial_positions.size())
    paramError("initial_radii", "Must have one entry per initial position");
  for (const auto r : _initial_radii)
    if (r <= 0)
      paramError("initial_radii", "Radii must be positive");
  for (const auto i : _initial_frozen)
    if (i >= _initial_positions.size())
      paramError("initial_frozen", "Index ", i, " beyond the initial positions");
  if (_insertion_radius_range.size() &&
      (_insertion_radius_range.size() != 2 || _insertion_radius_range[0] <= 0 ||
       _insertion_radius_range[1] < _insertion_radius_range[0] ||
       _insertion_radius_range[1] > _radius))
    paramError("insertion_radius_range",
               "Two radii, positive, in order, and no larger than radius (the neighbor list is "
               "sized by radius)");
  if (_initial_velocities.size() && _initial_velocities.size() != _initial_positions.size())
    paramError("initial_velocities", "Must have one entry per initial position");
  if (getParam<std::vector<Point>>("wall_points").size() !=
      getParam<std::vector<RealVectorValue>>("wall_normals").size())
    paramError("wall_normals", "Must have one entry per wall point");
  for (const auto & normal : getParam<std::vector<RealVectorValue>>("wall_normals"))
    if (normal.norm() == 0)
      paramError("wall_normals", "Normals must be nonzero");
  if ((_walls.n || _wall_boundaries.size()) && !contactEnabled())
    paramError("wall_points", "Walls need a contact model with a positive stiffness or modulus");
  for (const auto & normal : getParam<std::vector<RealVectorValue>>("wall_normals"))
    _wall_normals.push_back(normal.unit());
  if (_wall_velocities.empty())
    _wall_velocities.assign(_walls.n, RealVectorValue(0, 0, 0));
  else if (_wall_velocities.size() != _walls.n)
    paramError("wall_velocities", "One per wall point");
  if (_servo_walls.size() != _servo_forces.size())
    paramError("servo_forces", "One per servo wall");
  for (const auto w : _servo_walls)
    if (w >= _walls.n)
      paramError("servo_walls", "Wall ", w, " does not exist");
  if (_servo_walls.size() && (_servo_gain <= 0 || _servo_max_velocity <= 0))
    paramError("servo_gain", "Servo walls need a positive gain and maximum velocity");
  _walls.setVelocities(_wall_velocities);
  _walls.setRotation(getParam<Point>("wall_rotation_center"),
                     getParam<RealVectorValue>("wall_angular_velocity"));
  if (_insertion_box.size() && _insertion_box.size() != 2)
    paramError("insertion_box", "Two corners");
  if (_insertion_rate > 0 && _insertion_box.empty())
    paramError("insertion_rate", "Needs an insertion_box");
  if (_outflow_box.size() && _outflow_box.size() != 2)
    paramError("outflow_box", "Two corners");
  for (const auto & name : getParam<std::vector<VariableName>>("wall_displacements"))
    _wall_displacements.push_back(&_fe_problem.getStandardVariable(0, name));
  if (_wall_displacements.size() &&
      _wall_displacements.size() != _fe_problem.mesh().dimension())
    paramError("wall_displacements", "One variable per mesh dimension");
  if (_wall_displacements.size() && _wall_boundaries.empty())
    paramError("wall_displacements", "Only sideset walls (wall_boundaries) can move");
#ifndef DEM_HAVE_ARBORX
  if (_broad_phase == DEM::BroadPhase::ARBORX_BVH)
    paramError("broad_phase", "The ArborX broad phase needs the ArborX submodule (modules/dem/contrib/arborx)");
#endif
  // The far side of a periodic face is found by the point locator, which only sees the local and
  // ghost elements of a distributed mesh
  if (_periodic.size() && _fe_problem.mesh().isDistributedMesh())
    paramError("periodic", "Periodic directions require a replicated mesh");

  // Side centroids and normals are only cached on request
  _fe_problem.mesh().getKokkosMesh()->setNeedsElementSideGeometry();

  _checkpoint.cloud = this;
}

DEM::Hertz
KokkosParticleCloud::makeHertz() const
{
  if (_contact_model != ContactModel::HERTZ)
    return {};
  if (!isParamValid("youngs_modulus") || !isParamValid("poissons_ratio"))
    paramError("contact_model", "The Hertz model needs youngs_modulus and poissons_ratio");
  const Real e = getParam<Real>("restitution");
  const Real log_e = std::log(e);
  const Real E = getParam<Real>("youngs_modulus"), nu = getParam<Real>("poissons_ratio");
  // Two bodies of the same material: E* = E / (2 (1 - nu^2)) and G* = G / (2 (2 - nu)) with
  // G = E / (2 (1 + nu)); a given tangential_stiffness is 8 G* directly (LAMMPS's k_t)
  const Real k_t = getParam<Real>("tangential_stiffness");
  return {E / (2 * (1 - nu * nu)),
          -log_e / std::sqrt(log_e * log_e + libMesh::pi * libMesh::pi),
          k_t > 0 ? k_t / 8 : E / (2 * (1 + nu)) / (2 * (2 - nu)),
          getParam<Real>("dissipation_time"),
          getParam<bool>("limit_damping"),
          getParam<Real>("normal_damping"),
          getParam<Real>("tangential_damping")};
}

Real
KokkosParticleCloud::particleMass(const Real radius) const
{
  return _disk_mass ? _density * libMesh::pi * radius * radius
                    : _density * 4.0 / 3.0 * libMesh::pi * Utility::pow<3>(radius);
}

bool
KokkosParticleCloud::contactEnabled() const
{
  return _contact_model == ContactModel::HERTZ ? _hertz.enabled() : _linear.enabled();
}

void
KokkosParticleCloud::initialSetup()
{
  const auto & mesh = kokkosMesh();
  const auto num_total = mesh.getNumLocalAndPossiblyOneNeighborLayerGhostElements();

  // Reverse element ID map, owner ranks, and libMesh IDs for local and one-layer ghost elements
  _cid_to_elem.assign(num_total, nullptr);
  auto owner_host = ::Kokkos::View<processor_id_type *, ::Kokkos::HostSpace>("dem_owner", num_total);
  auto id_host = ::Kokkos::View<dof_id_type *, ::Kokkos::HostSpace>("dem_libmesh_id", num_total);
  for (const auto elem : _fe_problem.mesh().getMesh().active_local_element_ptr_range())
  {
    const auto cid = mesh.getContiguousElementID(elem);
    _cid_to_elem[cid] = elem;
    owner_host(cid) = processor_id();
    id_host(cid) = elem->id();
  }
  for (const auto & [elem, cid] : mesh.getGhostElemIdMapping())
  {
    _cid_to_elem[cid] = elem;
    owner_host(cid) = elem->processor_id();
    id_host(cid) = elem->id();
  }
  _owner = ::Kokkos::View<processor_id_type *>("dem_owner", num_total);
  _libmesh_id = ::Kokkos::View<dof_id_type *>("dem_libmesh_id", num_total);
  ::Kokkos::deep_copy(_owner, owner_host);
  ::Kokkos::deep_copy(_libmesh_id, id_host);

  _locator = _fe_problem.mesh().getPointLocator();
  _locator->enable_out_of_mesh_mode();

  if (!_checkpoint.loaded)
    placeInitialParticles();
  _next_gid = _initial_positions.size();

  // Ghosting is decided against the neighboring ranks' bounding boxes, gathered once since the
  // mesh does not move
  _r_max = _initial_radii.size() ? *std::max_element(_initial_radii.begin(), _initial_radii.end())
                                 : _radius;
  // A rank without elements has an invalid (inverted) local bounding box, which the inflation
  // turns into an infinite one that every rank would ghost to; keep it inverted so it neighbors
  // nothing
  const auto my_box = mesh.getNumLocalElements()
                          ? _fe_problem.mesh().getInflatedProcessorBoundingBox()
                          : libMesh::BoundingBox();
  std::vector<Real> bounds = {my_box.min()(0),
                              my_box.min()(1),
                              my_box.min()(2),
                              my_box.max()(0),
                              my_box.max()(1),
                              my_box.max()(2)};
  _communicator.allgather(bounds);
  _rank_boxes.clear();
  for (const auto rank : make_range(n_processors()))
    _rank_boxes.emplace_back(Point(bounds[6 * rank], bounds[6 * rank + 1], bounds[6 * rank + 2]),
                             Point(bounds[6 * rank + 3], bounds[6 * rank + 4], bounds[6 * rank + 5]));
  // The domain is the mesh bounding box; a periodic direction wraps over its extent
  const auto domain = libMesh::MeshTools::create_bounding_box(_fe_problem.mesh().getMesh());
  _domain_volume = 1;
  for (const auto c : make_range(3))
  {
    _domain_min(c) = domain.min()(c);
    _period(c) = _periodic.contains(std::string(1, "xyz"[c])) ? domain.max()(c) - domain.min()(c) : 0;
    if (c < _fe_problem.mesh().dimension())
      _domain_volume *= domain.max()(c) - domain.min()(c);
  }

  // Ghosts go to every rank, this one included, whose box our particles come within the cutoff
  // of under some combination of period shifts (-L, 0, +L per periodic direction), ordered by
  // rank so the ghosts sent to one rank are contiguous
  const Real cutoff = 2 * _r_max + _skin;
  _neighbor_ranks.clear();
  _ghost_targets.clear();

  // Sideset walls of the local and ghost elements; a face is a candidate of the elements within
  // a particle's reach of it
  if (_wall_boundaries.size())
  {
    // A name the mesh does not know would otherwise be given a fresh, empty ID
    const auto & boundary_info = _fe_problem.mesh().getMesh().get_boundary_info();
    const auto known = boundary_info.get_boundary_ids();
    for (const auto & name : _wall_boundaries)
    {
      const auto id = boundary_info.get_id_by_name(name);
      const bool numeric = !name.empty() && std::all_of(name.begin(), name.end(), ::isdigit);
      if (id == libMesh::BoundaryInfo::invalid_id &&
          !(numeric && known.count(static_cast<BoundaryID>(std::stoi(name)))))
        paramError("wall_boundaries", "The mesh has no boundary named '", name, "'");
    }
    _sideset_walls.build(_fe_problem.mesh(),
                         _cid_to_elem,
                         mesh.getNumLocalElements(),
                         _fe_problem.mesh().getBoundaryIDs(_wall_boundaries, true),
                         _r_max + _skin + _wall_reach,
                         _wall_curvature,
                         _rank_boxes,
                         _communicator);
    // Start from the initially displaced positions, at rest
    if (_wall_displacements.size())
      _sideset_walls.move(displacedWallVertices(), 0);
  }
  for (const auto rank : make_range(n_processors()))
    for (int sz = -1; sz <= 1; ++sz)
      for (int sy = -1; sy <= 1; ++sy)
        for (int sx = -1; sx <= 1; ++sx)
        {
          const Moose::Kokkos::Real3 shift(sx * _period(0), sy * _period(1), sz * _period(2));
          if ((sx && !_period(0)) || (sy && !_period(1)) || (sz && !_period(2)) ||
              (rank == processor_id() && !sx && !sy && !sz))
            continue;
          auto reach = my_box;
          reach.first += Point(shift(0) - cutoff, shift(1) - cutoff, shift(2) - cutoff);
          reach.second += Point(shift(0) + cutoff, shift(1) + cutoff, shift(2) + cutoff);
          if (!reach.intersects(_rank_boxes[rank]))
            continue;
          if (_neighbor_ranks.empty() || _neighbor_ranks.back() != rank)
            _neighbor_ranks.push_back(rank);
          _ghost_targets.push_back({_neighbor_ranks.size() - 1, shift});
        }

  checkMesh();
  checkTimestep();
  _set_up = true;

  if (_checkpoint.loaded)
    restore();
  else
  {
    updateNeighborList();
    computeForces(false);
    ::Kokkos::deep_copy(_cloud.f_out, _cloud.f);
    ::Kokkos::deep_copy(_cloud.tau_out, _cloud.tau);
    // Reduce so the initial values are global before the first execution
    count();
    finalize();
  }
}

void
KokkosParticleCloud::checkpoint(DEM::CloudCheckpoint & checkpoint)
{
  constexpr auto record_size = DEM::ParticleCloud::record_size;
  constexpr auto invalid_id = libMesh::DofObject::invalid_id;
  const auto cloud = _cloud;
  const auto n = _cloud.size();
  const auto num_elems = _libmesh_id.extent(0);
  const auto libmesh_id = _libmesh_id;

  // Particles as migration records, with the libMesh ID of the element, or none for a particle
  // that is unresolved (the sentinels are beyond the element count) so it is located on restore
  ::Kokkos::View<Real *> records("dem_checkpoint_particles", n * record_size);
  ::Kokkos::parallel_for(
      "dem_checkpoint_pack", n, KOKKOS_LAMBDA(const std::size_t i) {
        const auto elem = cloud.elem(i);
        cloud.pack(i, elem < num_elems ? libmesh_id(elem) : invalid_id, &records(i * record_size));
      });
  auto records_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, records);
  checkpoint.particles.assign(records_host.data(), records_host.data() + n * record_size);

  // Contact histories: count, scan, fill over the map's slots
  const auto states = _pair_states;
  const auto capacity = states.capacity();
  ::Kokkos::View<std::size_t *> slot("dem_checkpoint_slot", capacity + 1);
  std::size_t total = 0;
  ::Kokkos::parallel_scan(
      "dem_checkpoint_scan",
      capacity,
      KOKKOS_LAMBDA(const std::size_t e, std::size_t & sum, const bool final) {
        if (final)
          slot(e) = sum;
        sum += states.valid_at(e);
      },
      total);
  ::Kokkos::View<Real *> history_records("dem_checkpoint_histories",
                                          total * pair_state_record_size);
  ::Kokkos::parallel_for(
      "dem_checkpoint_histories", capacity, KOKKOS_LAMBDA(const std::size_t e) {
        if (!states.valid_at(e))
          return;
        const auto & key = states.key_at(e);
        const auto & value = states.value_at(e);
        Real * record = &history_records(slot(e) * pair_state_record_size);
        record[0] = key.lower;
        record[1] = key.upper;
        for (unsigned int c = 0; c < 3; ++c)
        {
          record[2 + c] = value.delta_t[c];
          record[5 + c] = value.delta_r[c];
        }
      });
  auto histories_host =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, history_records);
  checkpoint.histories.assign(histories_host.data(),
                              histories_host.data() + total * pair_state_record_size);

  checkpoint.analytic_wall_points = _walls.savePoints();
  if (_sideset_walls.n)
    _sideset_walls.save(checkpoint.wall_vertices, checkpoint.wall_velocities);
  else
  {
    checkpoint.wall_vertices.clear();
    checkpoint.wall_velocities.clear();
  }

  checkpoint.num_exited = _num_exited_local;
  checkpoint.num_migrated = _num_migrated;
  checkpoint.num_neighbor_list_builds = _num_neighbor_list_builds;
  checkpoint.num_ghost_forwards = _num_ghost_forwards;
  checkpoint.timestep_warned = _timestep_warned;
  checkpoint.insertion_draws = _insertion_draws;
  checkpoint.insertion_carry = _insertion_carry;
  checkpoint.next_gid = _next_gid;
}

void
KokkosParticleCloud::restore()
{
  constexpr auto record_size = DEM::ParticleCloud::record_size;
  constexpr auto invalid_id = libMesh::DofObject::invalid_id;
  const auto & mesh = kokkosMesh();
  const auto & libmesh_mesh = _fe_problem.mesh().getMesh();
  auto & checkpoint = _checkpoint;
  checkpoint.loaded = false;

  // Particles: each record's element must be one of this rank's, so the checkpoint must come
  // from the same mesh, partitioning, and number of ranks
  const auto n = checkpoint.particles.size() / record_size;
  DEM::ParticleCloud::HostMirror host(n);
  for (const auto j : make_range(n))
  {
    const Real * record = checkpoint.particles.data() + j * record_size;
    host.unpack(record, j);
    const auto id = static_cast<dof_id_type>(record[record_size - 1]);
    const auto * elem = id == invalid_id
                            ? (*_locator)(Point(host.x(j, 0), host.x(j, 1), host.x(j, 2)))
                            : libmesh_mesh.query_elem_ptr(id);
    if (!elem || elem->processor_id() != processor_id())
      mooseError("Restoring particle ",
                 host.gid(j),
                 " from the checkpoint failed: its element is not on this rank. The mesh, its "
                 "partitioning, and the number of ranks must be those of the checkpoint");
    host.elem(j) = mesh.getContiguousElementID(elem);
  }
  _cloud.allocate(n);
  _cloud.copyFrom(host, 0);
  _cloud.n = n;
  allocateScratch(n);
  _num_ghosts = 0;

  // Histories
  _pair_states = DEM::PairStateMap(2 * checkpoint.histories.size() / pair_state_record_size + 128);
  unpackPairStates(checkpoint.histories);

  // Moving walls resume where they were, with the velocities they had, so the forces below
  // see the same relative velocities as before the checkpoint
  if (checkpoint.analytic_wall_points.size() == 6 * _walls.n)
    _walls.loadPoints(checkpoint.analytic_wall_points);
  if (_sideset_walls.n)
  {
    if (checkpoint.wall_vertices.size() != 9 * _sideset_walls.n)
      mooseError("The checkpoint's sideset walls do not match this rank's");
    _sideset_walls.load(checkpoint.wall_vertices, checkpoint.wall_velocities);
  }
  checkpoint.wall_vertices.clear();
  checkpoint.wall_velocities.clear();

  _num_exited_local = checkpoint.num_exited;
  _num_migrated = checkpoint.num_migrated;
  _num_neighbor_list_builds = checkpoint.num_neighbor_list_builds;
  _num_ghost_forwards = checkpoint.num_ghost_forwards;
  _timestep_warned = checkpoint.timestep_warned;
  // The generator is replayed to where it was
  _insertion_generator.seed(getParam<unsigned int>("insertion_seed"));
  _insertion_generator.discard(checkpoint.insertion_draws);
  _insertion_draws = checkpoint.insertion_draws;
  _insertion_carry = checkpoint.insertion_carry;
  _next_gid = checkpoint.next_gid;
  checkpoint.particles.clear();
  checkpoint.histories.clear();

  // Rebuild everything derived: the neighbor list (and with it the ghosts and the pair-state
  // map, which keeps the restored histories of the listed pairs), the output forces, and the
  // counters; the restored forces are the last substep's and feed the next half-kick as they
  // would have without the restart
  _force_rebuild = true;
  updateNeighborList();
  std::swap(_cloud.f, _cloud.f_out);
  std::swap(_cloud.tau, _cloud.tau_out);
  computeForces(false);
  std::swap(_cloud.f, _cloud.f_out);
  std::swap(_cloud.tau, _cloud.tau_out);
  count();
  _initial_finalize = true;
  finalize();
}

void
DEM::dataStore(std::ostream & stream, DEM::CloudCheckpoint & checkpoint, void * context)
{
  mooseAssert(checkpoint.cloud, "Checkpoint without its cloud");
  checkpoint.cloud->checkpoint(checkpoint);
  ::dataStore(stream, checkpoint.particles, context);
  ::dataStore(stream, checkpoint.histories, context);
  ::dataStore(stream, checkpoint.wall_vertices, context);
  ::dataStore(stream, checkpoint.analytic_wall_points, context);
  ::dataStore(stream, checkpoint.wall_velocities, context);
  ::dataStore(stream, checkpoint.num_exited, context);
  ::dataStore(stream, checkpoint.num_migrated, context);
  ::dataStore(stream, checkpoint.num_neighbor_list_builds, context);
  ::dataStore(stream, checkpoint.num_ghost_forwards, context);
  ::dataStore(stream, checkpoint.timestep_warned, context);
  ::dataStore(stream, checkpoint.insertion_draws, context);
  ::dataStore(stream, checkpoint.insertion_carry, context);
  ::dataStore(stream, checkpoint.next_gid, context);
}

void
DEM::dataLoad(std::istream & stream, DEM::CloudCheckpoint & checkpoint, void * context)
{
  ::dataLoad(stream, checkpoint.particles, context);
  ::dataLoad(stream, checkpoint.histories, context);
  ::dataLoad(stream, checkpoint.wall_vertices, context);
  ::dataLoad(stream, checkpoint.analytic_wall_points, context);
  ::dataLoad(stream, checkpoint.wall_velocities, context);
  ::dataLoad(stream, checkpoint.num_exited, context);
  ::dataLoad(stream, checkpoint.num_migrated, context);
  ::dataLoad(stream, checkpoint.num_neighbor_list_builds, context);
  ::dataLoad(stream, checkpoint.num_ghost_forwards, context);
  ::dataLoad(stream, checkpoint.timestep_warned, context);
  ::dataLoad(stream, checkpoint.insertion_draws, context);
  ::dataLoad(stream, checkpoint.insertion_carry, context);
  ::dataLoad(stream, checkpoint.next_gid, context);
  checkpoint.loaded = true;
  // After setup (a restep or backup) restore at once; before it, initialSetup() will
  if (checkpoint.cloud && checkpoint.cloud->isSetUp())
    checkpoint.cloud->restore();
}

void
KokkosParticleCloud::placeInitialParticles()
{
  // Each rank keeps the particles that start in one of its local elements
  const auto & mesh = kokkosMesh();
  std::vector<std::pair<int64_t, ContiguousElementID>> local;
  for (const auto i : index_range(_initial_positions))
  {
    const auto * elem = (*_locator)(_initial_positions[i]);
    if (elem && elem->processor_id() == processor_id())
      local.emplace_back(i, mesh.getContiguousElementID(elem));
  }

  const auto n = local.size();
  DEM::ParticleCloud::HostMirror host(n);
  for (const auto j : make_range(n))
  {
    const auto [gid, cid] = local[j];
    const Real radius = _initial_radii.size() ? _initial_radii[gid] : _radius;
    const Real mass = particleMass(radius);
    // Solid sphere
    const Real inertia = 0.4 * mass * radius * radius;
    host.gid(j) = gid;
    for (const auto c : make_range(3))
    {
      host.x(j, c) = _initial_positions[gid](c);
      host.v(j, c) = _initial_velocities.size() ? _initial_velocities[gid](c) : _initial_velocity(c);
      host.omega(j, c) = _initial_angular_velocity(c);
      host.f(j, c) = 0;
      host.tau(j, c) = 0;
    }
    // Identity orientation
    host.q(j, 0) = 1;
    host.q(j, 1) = host.q(j, 2) = host.q(j, 3) = 0;
    host.r(j) = radius;
    host.m(j) = mass;
    host.inertia(j) = inertia;
    host.elem(j) = cid;
    host.target_rank(j) = libMesh::DofObject::invalid_processor_id;
    host.hops(j) = 0;
    host.frozen(j) = _initial_frozen.count(gid) ? 1 : 0;
  }
  _cloud.allocate(n);
  _cloud.copyFrom(host, 0);
  _cloud.n = n;
  allocateScratch(n);
}

void
KokkosParticleCloud::servoWalls()
{
  // Proportional control on the force measured at the end of the last step: the compressive
  // force is the particles' push against the wall's inward normal. The velocity is capped so
  // the wall cannot change the force of a single contact by more than the error within the
  // step, with the contact stiffness at the target force (the servo sees the force one step
  // late), which keeps the control stable for any gain
  const Real dt = _dt;
  for (const auto k : index_range(_servo_walls))
  {
    const auto w = _servo_walls[k];
    const auto normal = _wall_normals[w];
    const Real target = _servo_forces[k];
    const Real error = target - (-(wallForce(w) * normal));
    Real stiffness = _linear.stiffness;
    if (_contact_model == ContactModel::HERTZ)
    {
      const Real r = _r_max;
      const Real overlap = std::pow(std::max(target, 1e-12) / (4.0 / 3.0 * _hertz.effective_modulus * std::sqrt(r)), 2.0 / 3.0);
      stiffness = 2 * _hertz.effective_modulus * std::sqrt(r * overlap);
    }
    Real limit = _servo_max_velocity;
    if (stiffness > 0 && dt > 0)
      limit = std::min(limit, std::abs(error) / (stiffness * dt));
    const Real v = std::clamp(_servo_gain * error, -limit, limit);
    _wall_velocities[w] = v * normal;
  }
  _walls.setVelocities(_wall_velocities);
}

void
KokkosParticleCloud::appendParticles(const DEM::ParticleCloud::HostMirror & host,
                                     const std::size_t count)
{
  const auto begin = _cloud.size();
  const auto end = begin + count;
  if (end > _cloud.capacity())
    reserve(std::max(end, 2 * _cloud.capacity()));
  _cloud.copyFrom(host, begin);
  _cloud.n = end;
}

void
KokkosParticleCloud::insertParticles()
{
  if (_t < _insertion_start_time || _t > _insertion_end_time)
    return;
  const auto & mesh = kokkosMesh();
  const Point lower(std::min(_insertion_box[0](0), _insertion_box[1](0)),
                    std::min(_insertion_box[0](1), _insertion_box[1](1)),
                    std::min(_insertion_box[0](2), _insertion_box[1](2)));
  const Point upper(std::max(_insertion_box[0](0), _insertion_box[1](0)),
                    std::max(_insertion_box[0](1), _insertion_box[1](1)),
                    std::max(_insertion_box[0](2), _insertion_box[1](2)));
  const Real radius = _radius;

  // Whole particles this step, the fraction carried over
  _insertion_carry += _insertion_rate * _dt;
  const auto wanted = static_cast<std::size_t>(std::floor(_insertion_carry));
  if (!wanted)
    return;

  // Positions of the local particles near the box, on host, to test candidates against
  const auto cloud = _cloud;
  const Moose::Kokkos::Real3 lo(lower(0) - 2 * _r_max, lower(1) - 2 * _r_max, lower(2) - 2 * _r_max);
  const Moose::Kokkos::Real3 hi(upper(0) + 2 * _r_max, upper(1) + 2 * _r_max, upper(2) + 2 * _r_max);
  ::Kokkos::View<std::size_t *> slot("dem_insert_slot", cloud.size() + 1);
  std::size_t near = 0;
  ::Kokkos::parallel_scan(
      "dem_insert_near",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum, const bool final) {
        bool inside = cloud.elem(i) != DEM::ParticleCloud::exited;
        for (unsigned int c = 0; c < 3; ++c)
          inside = inside && cloud.x(i, c) >= lo(c) && cloud.x(i, c) <= hi(c);
        if (final)
          slot(i) = inside ? sum : std::numeric_limits<std::size_t>::max();
        sum += inside;
      },
      near);
  ::Kokkos::View<Real * [4], ::Kokkos::LayoutRight> packed("dem_insert_near_particles", near);
  ::Kokkos::parallel_for(
      "dem_insert_pack", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        const auto s = slot(i);
        if (s == std::numeric_limits<std::size_t>::max())
          return;
        for (unsigned int c = 0; c < 3; ++c)
          packed(s, c) = cloud.x(i, c);
        packed(s, 3) = cloud.r(i);
      });
  auto existing = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, packed);

  // Every rank draws the same candidates; a candidate is rejected when it overlaps a particle
  // on any rank (a global decision) or an earlier accepted candidate, with the mesh dimension's
  // coordinates only, and at most ten tries per particle wanted
  // One engine call per coordinate, so discard() replays the sequence exactly on restart
  const auto unit = [this]()
  {
    ++_insertion_draws;
    return (_insertion_generator() >> 11) * 0x1.0p-53;
  };
  const auto dim = _fe_problem.mesh().dimension();
  std::vector<Point> accepted;
  std::vector<Real> accepted_radii;
  for (unsigned int attempt = 0; attempt < 10 * wanted && accepted.size() < wanted; ++attempt)
  {
    Point p;
    for (const auto c : make_range(3))
    {
      const Real u = unit();
      p(c) = c < dim ? lower(c) + u * (upper(c) - lower(c)) : lower(c);
    }
    const Real rp = _insertion_radius_range.size()
                        ? _insertion_radius_range[0] +
                              unit() * (_insertion_radius_range[1] - _insertion_radius_range[0])
                        : radius;
    unsigned int overlaps = 0;
    for (const auto k : index_range(accepted))
      overlaps = overlaps || (p - accepted[k]).norm() < rp + accepted_radii[k];
    for (const auto e : make_range(near))
      overlaps = overlaps ||
                 (p - Point(existing(e, 0), existing(e, 1), existing(e, 2))).norm() <
                     rp + existing(e, 3);
    _communicator.max(overlaps);
    if (!overlaps)
    {
      accepted.push_back(p);
      accepted_radii.push_back(rp);
    }
  }
  _insertion_carry -= accepted.size();

  // Each rank keeps the accepted particles that lie in its elements
  std::vector<std::pair<Point, ContiguousElementID>> mine;
  std::vector<int64_t> gids;
  std::vector<Real> radii;
  for (const auto k : index_range(accepted))
  {
    const auto * elem = (*_locator)(accepted[k]);
    if (elem && elem->processor_id() == processor_id())
    {
      mine.emplace_back(accepted[k], mesh.getContiguousElementID(elem));
      gids.push_back(_next_gid);
      radii.push_back(accepted_radii[k]);
    }
    ++_next_gid;
  }
  if (mine.empty())
    return;
  DEM::ParticleCloud::HostMirror host(mine.size());
  for (const auto j : index_range(mine))
  {
    const Real radius = radii[j];
    const Real mass = particleMass(radius);
    host.gid(j) = gids[j];
    for (const auto c : make_range(3))
    {
      host.x(j, c) = mine[j].first(c);
      host.v(j, c) = _insertion_velocity(c);
      host.omega(j, c) = 0;
      host.f(j, c) = 0;
      host.tau(j, c) = 0;
    }
    host.q(j, 0) = 1;
    host.q(j, 1) = host.q(j, 2) = host.q(j, 3) = 0;
    host.r(j) = radius;
    host.m(j) = mass;
    host.inertia(j) = 0.4 * mass * radius * radius;
    host.elem(j) = mine[j].second;
    host.target_rank(j) = libMesh::DofObject::invalid_processor_id;
    host.hops(j) = 0;
    host.frozen(j) = 0;
  }
  appendParticles(host, mine.size());
}

void
KokkosParticleCloud::applyOutflow()
{
  const auto cloud = _cloud;
  const Moose::Kokkos::Real3 lo(std::min(_outflow_box[0](0), _outflow_box[1](0)),
                                std::min(_outflow_box[0](1), _outflow_box[1](1)),
                                std::min(_outflow_box[0](2), _outflow_box[1](2)));
  const Moose::Kokkos::Real3 hi(std::max(_outflow_box[0](0), _outflow_box[1](0)),
                                std::max(_outflow_box[0](1), _outflow_box[1](1)),
                                std::max(_outflow_box[0](2), _outflow_box[1](2)));
  constexpr auto exited = DEM::ParticleCloud::exited;
  ::Kokkos::parallel_for(
      "dem_outflow", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        bool inside = true;
        for (unsigned int c = 0; c < 3; ++c)
          inside = inside && cloud.x(i, c) >= lo(c) && cloud.x(i, c) <= hi(c);
        if (inside)
          cloud.elem(i) = exited;
      });
}

void
KokkosParticleCloud::checkMesh()
{
  // The face walk assumes planar faces: every side's vertices must lie in the plane of its
  // first three, to a tolerance relative to the element size
  Real worst = 0;
  Real h_min = std::numeric_limits<Real>::max();
  for (const auto * elem : _cid_to_elem)
  {
    const Real h = elem->hmin();
    h_min = std::min(h_min, h);
    for (const auto s : make_range(elem->n_sides()))
    {
      const auto side = elem->side_ptr(s);
      if (side->n_vertices() < 4)
        continue;
      const Point a = side->point(0), b = side->point(1), c = side->point(2);
      Point normal = (b - a).cross(c - a);
      if (normal.norm() == 0)
        continue;
      normal /= normal.norm();
      for (const auto k : make_range(3u, side->n_vertices()))
        worst = std::max(worst, std::abs((side->point(k) - a) * normal) / h);
    }
  }
  _communicator.max(worst);
  _communicator.min(h_min);
  // 1e-8: roundoff of rotated or transformed meshes is 1e-16, a warped element far larger
  if (worst > 1e-8)
    mooseError("The mesh has a face that is not planar (a vertex ",
               worst,
               " element sizes out of its plane); the face walk needs planar faces");

  libmesh_ignore(h_min);
}

void
KokkosParticleCloud::checkTimestep()
{
  if (_timestep_check == "none" || !contactEnabled())
    return;
  const Real dt = _dt / _substeps;
  const Real r_min = _initial_radii.size()
                         ? *std::min_element(_initial_radii.begin(), _initial_radii.end())
                         : _radius;
  const Real m_min = particleMass(r_min);
  // A pair of the smallest particles has the shortest contact time
  const Real m_eff = m_min / 2, r_eff = r_min / 2;
  std::vector<std::pair<std::string, Real>> scales;
  if (_contact_model == ContactModel::HERTZ)
  {
    // Rayleigh time of the smallest particle (LIGGGHTS check/timestep/gran)
    const Real E = getParam<Real>("youngs_modulus"), nu = getParam<Real>("poissons_ratio");
    const Real G = E / (2 * (1 + nu));
    scales.emplace_back("Rayleigh",
                        libMesh::pi * r_min * std::sqrt(_density / G) / (0.1631 * nu + 0.8766));
    // Hertz contact time at the largest relative speed, twice the largest particle speed
    const auto cloud = _cloud;
    Real v_max_sq = 0;
    ::Kokkos::parallel_reduce(
        "dem_check_vmax",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & max) {
          const Real v_sq = libMesh::Kokkos::make_vector_ref(cloud.v, i).norm_sq();
          max = v_sq > max ? v_sq : max;
        },
        ::Kokkos::Max<Real>(v_max_sq));
    _communicator.max(v_max_sq);
    if (v_max_sq > 0)
      scales.emplace_back(
          "Hertz",
          2.87 * std::pow(m_eff * m_eff /
                              (r_eff * Utility::pow<2>(_hertz.effective_modulus) * 2 *
                               std::sqrt(v_max_sq)),
                          0.2));
  }
  else
  {
    // Contact duration of the linear spring: half a period of the stiffest spring
    const Real k = std::max(_linear.stiffness, _linear.tangential_stiffness);
    scales.emplace_back("linear-spring contact", libMesh::pi * std::sqrt(m_eff / k));
  }

  for (const auto & [name, scale] : scales)
    if (dt > _timestep_fraction * scale)
    {
      const auto message = "The substep " + Moose::stringify(dt) + " is " +
                           Moose::stringify(100 * dt / scale) + "% of the " + name + " time " +
                           Moose::stringify(scale) + ", above the allowed " +
                           Moose::stringify(100 * _timestep_fraction) + "%";
      if (_timestep_check == "error")
        mooseError(message);
      else if (!_timestep_warned)
      {
        mooseWarning(message, "; increase substeps (further violations are not reported)");
        _timestep_warned = true;
      }
    }
}

void
KokkosParticleCloud::meshChanged()
{
  mooseError("The mesh changed; the particle cloud does not support adaptivity or mesh "
             "modification after setup");
}

std::vector<Point>
KokkosParticleCloud::displacedWallVertices() const
{
  const auto & mesh = _fe_problem.mesh().getMesh();
  const auto num = 3 * _sideset_walls.n;
  std::vector<Point> positions;
  positions.reserve(num);
  for (const auto & nodes : _sideset_walls.nodes)
    for (const auto id : nodes)
      positions.push_back(mesh.node_ref(id));

  // The nodes of ghost elements need not have their displacement dofs ghosted, so the values
  // are gathered collectively by index (every rank calls this once per step)
  for (const auto c : index_range(_wall_displacements))
  {
    const auto & var = *_wall_displacements[c];
    std::vector<libMesh::dof_id_type> dofs;
    dofs.reserve(num);
    for (const auto & nodes : _sideset_walls.nodes)
      for (const auto id : nodes)
        dofs.push_back(mesh.node_ref(id).dof_number(var.sys().number(), var.number(), 0));
    std::vector<Number> values;
    var.sys().currentSolution()->localize(values, dofs);
    for (const auto k : make_range(num))
      positions[k](c) += values[k];
  }
  return positions;
}

void
KokkosParticleCloud::compute()
{
  const auto cloud = _cloud;
  ::Kokkos::parallel_for(
      "dem_reset_hops", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) { cloud.hops(i) = 0; });

  // Moving walls travel to this step's displaced positions over the step
  if (_wall_displacements.size())
    _sideset_walls.move(displacedWallVertices(), _dt);
  if (_insertion_rate > 0)
    insertParticles();
  if (_servo_walls.size())
    servoWalls();
  // The Hertz time depends on the current speeds
  checkTimestep();

  // Velocity Verlet; the forces at the start of the step were computed at the end of the last one
  _num_migrated = 0;
  for (unsigned int substep = 0; substep < _substeps; ++substep)
  {
    kickDrift();
    _sideset_walls.advance(_dt / _substeps);
    _walls.advance(_dt / _substeps);
    walk();
    // A particle that has strayed beyond its rank's box, across a partition or a periodic face,
    // is settled at once: beyond the skin its far-side neighbors are no longer ghosted here
    if (strayed())
      settle();
    if (!updateNeighborList())
      startGhostUpdate();
    computeForces(true);
    kick();
  }

  // Every particle is settled at the end of the step, so it is located and owned where it is;
  // those in the outflow box leave with the exited. The forces of the last substep travel with
  // the particles for the next step's first half-kick
  if (_outflow_box.size())
    applyOutflow();
  settle();
  // The forces for output, contact counts, and wall forces are recomputed from the settled
  // positions and the full-step velocities into the output views, leaving the half-kick's
  // forces (computed with the half-step velocities, as in every substep) untouched; the
  // histories were advanced for this substep already and travelled with the particles
  if (!updateNeighborList())
    startGhostUpdate();
  std::swap(_cloud.f, _cloud.f_out);
  std::swap(_cloud.tau, _cloud.tau_out);
  computeForces(false);
  std::swap(_cloud.f, _cloud.f_out);
  std::swap(_cloud.tau, _cloud.tau_out);

  count();
  // Every rank takes part in the reduction, whether or not it holds wall faces
  if (_wall_boundaries.size())
  {
    int overflow = 0;
    if (_sideset_walls.n)
      ::Kokkos::deep_copy(overflow, _sideset_walls.overflow);
    _communicator.max(overflow);
    if (overflow)
      mooseError("A particle had more than ",
                 DEM::SidesetWalls::max_contacts,
                 " candidate sideset wall contacts; the walls are too finely meshed for the "
                 "particle size");
  }
  if (_verify)
    verify();
  if (_verify_pair_states)
    verifyPairStates();
}

bool
KokkosParticleCloud::strayed()
{
  // Beyond skin / 2 outside the box a particle may, at the next rebuild, be within the cutoff of
  // a particle that is not ghosted here (ghosts come within the cutoff of the box, not of the
  // particles); the flag is global so every rank settles together
  const auto cloud = _cloud;
  const auto & box = _rank_boxes[processor_id()];
  const Moose::Kokkos::Real3 lower(box.min()(0) - 0.5 * _skin, box.min()(1) - 0.5 * _skin, box.min()(2) - 0.5 * _skin);
  const Moose::Kokkos::Real3 upper(box.max()(0) + 0.5 * _skin, box.max()(1) + 0.5 * _skin, box.max()(2) + 0.5 * _skin);
  constexpr auto exited = DEM::ParticleCloud::exited;
  unsigned int flag = 0;
  ::Kokkos::parallel_reduce(
      "dem_strayed",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, unsigned int & max) {
        if (cloud.elem(i) == exited)
          return;
        for (unsigned int c = 0; c < 3; ++c)
          if (cloud.x(i, c) < lower(c) || cloud.x(i, c) > upper(c))
            max = 1;
      },
      ::Kokkos::Max<unsigned int>(flag));
  if (n_processors() > 1)
    _communicator.max(flag);
  return flag;
}

void
KokkosParticleCloud::settle()
{
  // Particles the walk could not resolve get a second chance from the point locator, and may
  // turn out to belong to another rank; the walk that follows flags those for migration
  resolveUnresolved();
  walk();

  // Remove the exited particles and hand the ones that reached a ghost element to their owner,
  // who continues the walk (with the same fallback, since the hop budget may not reach from the
  // element the sender saw), until every particle has settled. Each round moves a particle at
  // least one partition, or through a periodic face, so the number of rounds is bounded; the
  // bound is generous because the walk, which exits through the face the point is furthest in
  // front of, can wander across a partition and back on an irregular mesh
  const unsigned int max_rounds = 4 * n_processors() * (1 + _periodic.size());
  for (unsigned int round = 0;; ++round)
  {
    auto n_sent = migrate();
    _num_migrated += n_sent;
    _communicator.sum(n_sent);
    if (!n_sent)
      break;
    if (round == max_rounds)
      mooseError("Particles were still crossing partitions after ", max_rounds, " rounds");
    walk();
    resolveUnresolved();
    walk();
  }
}

template <typename Model>
DEM::ContactEvaluator<Model>
KokkosParticleCloud::evaluator(const Model & contact) const
{
  return {_cloud,
          _walls,
          _sideset_walls,
          contact,
          _friction,
          _wall_friction,
          _pair_states,
          _dt / _substeps,
          _deformed_torque_arm};
}

template <typename Model>
void
KokkosParticleCloud::computeForces(const Model & contact,
                                   const ::Kokkos::View<std::size_t *> & subset,
                                   const std::size_t count)
{
  const auto cloud = _cloud;
  const Moose::Kokkos::Real3 g(_gravity(0), _gravity(1), _gravity(2));
  const bool enabled = contact.enabled();
  const auto eval = evaluator(contact);
  const auto num_walls = _walls.n;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  constexpr auto exited = DEM::ParticleCloud::exited;

  ::Kokkos::parallel_for(
      "dem_forces", count, KOKKOS_LAMBDA(const std::size_t s) {
        const auto i = subset(s);
        Moose::Kokkos::Real3 f = 0, tau = 0;
        // Exited particles are inert until they are removed
        if (cloud.elem(i) != exited)
        {
          f = cloud.m(i) * g;

          // Sphere-sphere narrow phase over the neighbor list; every pair is listed from both
          // sides so each particle sums its own contacts, reading the history advanced by
          // updatePairStates(). The tangential force acts at the contact point, -r_i n from the
          // center
          if (enabled)
          {
            DEM::Contact c;
            Real f_n;
            Moose::Kokkos::Real3 f_t, tau_r;
            for (auto k = offsets(i); k < offsets(i + 1); ++k)
            {
              const auto j = pairs(k);
              if (cloud.elem(j) == exited || !eval.sphere(i, j, c))
                continue;
              const Real sign = cloud.gid(i) < cloud.gid(j) ? 1 : -1;
              eval.forces(
                  c, DEM::PairKey::pair(cloud.gid(i), cloud.gid(j)), sign, false, f_n, f_t, tau_r);
              f += f_n * c.normal + f_t;
              tau += tau_r - c.arm * c.normal.cross_product(f_t);
            }

            // Sphere-plane narrow phase against every wall
            for (std::size_t w = 0; w < num_walls; ++w)
            {
              if (!eval.wall(i, w, c))
                continue;
              eval.forces(c, DEM::PairKey::wall(cloud.gid(i), w), 1, false, f_n, f_t, tau_r);
              f += f_n * c.normal + f_t;
              tau += tau_r - c.arm * c.normal.cross_product(f_t);
            }

            // Sphere-face narrow phase against the sideset walls near the particle's element
            DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
            const auto n_wc = eval.sidesetContacts(i, cloud.r(i), wc);
            for (unsigned int k = 0; k < n_wc; ++k)
            {
              c = eval.sideset(i, wc[k]);
              eval.forces(c, eval.sidesetKey(i, wc[k]), 1, false, f_n, f_t, tau_r);
              f += f_n * c.normal + f_t;
              tau += tau_r - c.arm * c.normal.cross_product(f_t);
            }
          }
        }

        for (unsigned int c = 0; c < 3; ++c)
        {
          cloud.f(i, c) = f(c);
          cloud.tau(i, c) = tau(c);
        }
      });
}

template <typename Model>
void
KokkosParticleCloud::updatePairStates(const Model & contact, const bool ghost_pairs)
{
  if (!contact.enabled() || !historyNeeded())
    return;
  const auto cloud = _cloud;
  const auto n_local = _cloud.size();
  const auto eval = evaluator(contact);
  const auto num_walls = _walls.n;
  const Real skin = _skin;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  const auto ghost_self = _ghost_self;
  constexpr auto exited = DEM::ParticleCloud::exited;

  // Each pair is advanced from exactly one side per rank: a pair of local particles from the
  // lower index, a pair with another rank's ghost from the local particle (the other rank
  // advances its own copy identically), and a pair with a periodic image of a local particle,
  // which is listed from both sides, from the lower global ID. A listed pair that does not
  // overlap has its history reset
  ::Kokkos::parallel_for(
      "dem_pair_states", n_local, KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited)
          return;
        DEM::Contact c;
        Real f_n;
        Moose::Kokkos::Real3 f_t, tau_r;
        for (auto k = offsets(i); k < offsets(i + 1); ++k)
        {
          const auto j = pairs(k);
          const bool ghost = j >= n_local;
          if (ghost != ghost_pairs || cloud.elem(j) == exited)
            continue;
          const bool lower = cloud.gid(i) < cloud.gid(j);
          if (ghost ? (ghost_self(j - n_local) && !lower) : i > j)
            continue;
          const auto key = DEM::PairKey::pair(cloud.gid(i), cloud.gid(j));
          if (eval.sphere(i, j, c))
            eval.forces(c, key, lower ? 1 : -1, true, f_n, f_t, tau_r);
          else
            eval.reset(key);
        }
        if (!ghost_pairs)
        {
          for (std::size_t w = 0; w < num_walls; ++w)
          {
            const auto key = DEM::PairKey::wall(cloud.gid(i), w);
            if (eval.wall(i, w, c))
              eval.forces(c, key, 1, true, f_n, f_t, tau_r);
            else
              eval.reset(key);
          }
          // The touching faces (the reduced contacts) inherit the history of a same-surface
          // neighbor the contact slid off, then advance; every other listed face (within the
          // skin, unreduced) is reset
          DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
          DEM::WallContact listed[DEM::SidesetWalls::max_contacts];
          const auto n_touching = eval.sidesetContacts(i, cloud.r(i), wc);
          const auto n_listed = eval.sidesetContacts(i, cloud.r(i) + skin, listed, false);
          for (unsigned int k = 0; k < n_touching; ++k)
            eval.transferHistory(i, wc, n_touching, k);
          for (unsigned int k = 0; k < n_touching; ++k)
          {
            c = eval.sideset(i, wc[k]);
            eval.forces(c, eval.sidesetKey(i, wc[k]), 1, true, f_n, f_t, tau_r);
          }
          for (unsigned int k = 0; k < n_listed; ++k)
          {
            bool touching = false;
            for (unsigned int l = 0; l < n_touching && !touching; ++l)
              touching = wc[l].face == listed[k].face;
            if (!touching)
              eval.reset(eval.sidesetKey(i, listed[k]));
          }
        }
      });
}

void
KokkosParticleCloud::updatePairStates(const bool ghost_pairs)
{
  if (_contact_model == ContactModel::HERTZ)
    updatePairStates(_hertz, ghost_pairs);
  else
    updatePairStates(_linear, ghost_pairs);
}

void
KokkosParticleCloud::computeForces(const ::Kokkos::View<std::size_t *> & subset,
                                   const std::size_t count)
{
  // The model is chosen once here so the kernel inlines it (plan decision D5)
  if (_contact_model == ContactModel::HERTZ)
    computeForces(_hertz, subset, count);
  else
    computeForces(_linear, subset, count);
}

void
KokkosParticleCloud::computeForces(const bool advance_histories)
{
  // Interior particles only have pairs of local particles, so their histories and forces do not
  // need the ghost update and overlap its communication
  if (advance_histories)
    updatePairStates(false);
  computeForces(_neighbor_list.interior(), _neighbor_list.numInterior());
  finishGhostUpdate();
  if (advance_histories)
    updatePairStates(true);
  computeForces(_neighbor_list.boundary(), _neighbor_list.numBoundary());
}

void
KokkosParticleCloud::kickDrift()
{
  const auto cloud = _cloud;
  const Real dt = _dt / _substeps;
  constexpr auto exited = DEM::ParticleCloud::exited;

  ::Kokkos::parallel_for(
      "dem_kick_drift", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited || cloud.frozen(i))
          return;
        const Real half_dt_over_m = 0.5 * dt / cloud.m(i);
        const Real half_dt_over_I = 0.5 * dt / cloud.inertia(i);
        auto x = libMesh::Kokkos::make_vector_ref(cloud.x, i);
        auto v = libMesh::Kokkos::make_vector_ref(cloud.v, i);
        auto omega = libMesh::Kokkos::make_vector_ref(cloud.omega, i);

        v.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.f, i), half_dt_over_m);
        x.add_scaled(v, dt);
        omega.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.tau, i), half_dt_over_I);

        // Rotate the orientation by the body-frame angular velocity over the substep, exactly for
        // constant angular velocity: q <- q * (cos(theta/2), sin(theta/2) * omega_hat)
        const Real omega_norm = omega.norm();
        if (omega_norm > 0)
        {
          const Real half_theta = 0.5 * omega_norm * dt;
          const Real s = ::Kokkos::sin(half_theta) / omega_norm;
          const Real dw = ::Kokkos::cos(half_theta);
          const Real dx = s * omega(0), dy = s * omega(1), dz = s * omega(2);
          const Real qw = cloud.q(i, 0), qx = cloud.q(i, 1), qy = cloud.q(i, 2),
                     qz = cloud.q(i, 3);
          cloud.q(i, 0) = qw * dw - qx * dx - qy * dy - qz * dz;
          cloud.q(i, 1) = qw * dx + qx * dw + qy * dz - qz * dy;
          cloud.q(i, 2) = qw * dy - qx * dz + qy * dw + qz * dx;
          cloud.q(i, 3) = qw * dz + qx * dy - qy * dx + qz * dw;
          // Renormalize to keep roundoff from accumulating in the unit quaternion
          const Real inv_norm =
              1.0 / ::Kokkos::sqrt(cloud.q(i, 0) * cloud.q(i, 0) + cloud.q(i, 1) * cloud.q(i, 1) +
                                   cloud.q(i, 2) * cloud.q(i, 2) + cloud.q(i, 3) * cloud.q(i, 3));
          for (unsigned int c = 0; c < 4; ++c)
            cloud.q(i, c) *= inv_norm;
        }
      });
}

void
KokkosParticleCloud::kick()
{
  const auto cloud = _cloud;
  const Real dt = _dt / _substeps;
  constexpr auto exited = DEM::ParticleCloud::exited;

  ::Kokkos::parallel_for(
      "dem_kick", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited || cloud.frozen(i))
          return;
        auto v = libMesh::Kokkos::make_vector_ref(cloud.v, i);
        auto omega = libMesh::Kokkos::make_vector_ref(cloud.omega, i);
        v.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.f, i), 0.5 * dt / cloud.m(i));
        omega.add_scaled(libMesh::Kokkos::make_vector_ref(cloud.tau, i),
                         0.5 * dt / cloud.inertia(i));
      });
}

void
KokkosParticleCloud::walk()
{
  // Shallow copies captured by value into the kernel
  const Moose::Kokkos::Mesh mesh = kokkosMesh();
  const auto num_local = mesh.getNumLocalElements();
  const auto cloud = _cloud;
  const auto owner = _owner;
  const auto max_hops = _max_hops;
  const bool periodic = _periodic.size();
  constexpr auto invalid_elem = libMesh::DofObject::invalid_id;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;
  constexpr auto exited = DEM::ParticleCloud::exited;
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;

  ::Kokkos::parallel_for(
      "dem_walk", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        auto elem = cloud.elem(i);
        // An unresolved particle may carry the target the point locator gave it
        if (elem == unresolved)
          return;
        cloud.target_rank(i) = invalid_proc;
        if (elem == exited)
          return;

        const Moose::Kokkos::Real3 xp(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
        unsigned int hops = 0;
        while (true)
        {
          // Ghost elements have no cached side geometry; the owner continues the walk after
          // migration
          if (elem >= num_local)
          {
            cloud.target_rank(i) = owner(elem);
            break;
          }

          // The particle is inside if it is behind every outward face plane; otherwise leave
          // through the face it is furthest in front of
          const auto & info = mesh.getElementInfo(elem);
          Real d_max = 0;
          unsigned int exit_side = libMesh::invalid_uint;
          // Raw loop: make_range() is not device-callable on GPU backends
          for (unsigned int side = 0; side < mesh.getNumSides(info.type); ++side)
          {
            const Real d =
                (xp - mesh.getSideCentroid(elem, side)).dot_product(mesh.getSideNormal(elem, side));
            if (d > d_max)
            {
              d_max = d;
              exit_side = side;
            }
          }
          if (exit_side == libMesh::invalid_uint)
            break;

          if (hops == max_hops)
          {
            elem = unresolved;
            break;
          }
          elem = mesh.getNeighbor(elem, exit_side);
          ++hops;
          if (elem == invalid_elem)
          {
            // Left through a face without a neighbor, i.e. the mesh boundary. With periodic
            // directions the particle may have gone through a periodic face: it keeps its
            // unwrapped position, which its neighbors and their ghost images see, until the end
            // of the step, when it is wrapped and placed by the point locator
            elem = periodic ? unresolved : exited;
            break;
          }
        }
        cloud.elem(i) = elem;
        cloud.hops(i) += hops;
      });
}

void
KokkosParticleCloud::resolveUnresolved()
{
  const auto cloud = _cloud;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;
  std::size_t n_unresolved = 0;
  ::Kokkos::parallel_reduce(
      "dem_count_unresolved",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        sum += cloud.elem(i) == unresolved;
      },
      n_unresolved);
  if (!n_unresolved)
    return;

  // Only the unresolved particles cross to the host; a particle found in a ghost element is
  // handed off by the walk that follows, one found beyond the ghost layer (across a periodic
  // face) stays unresolved but is sent to the owner by migrate() with its element's libMesh ID
  const auto & mesh = kokkosMesh();
  const auto & ghost_ids = mesh.getGhostElemIdMapping();
  auto x_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto elem_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.elem);
  auto target_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.target_rank);
  auto far_host = ::Kokkos::create_mirror_view(_far_elem);
  for (const auto i : make_range(_cloud.size()))
  {
    if (elem_host(i) != unresolved)
      continue;
    // Wrap through the periodic faces; a particle still outside the mesh then has exited
    for (const auto c : make_range(3))
      if (_period(c) > 0)
        x_host(i, c) = _domain_min(c) + std::fmod(x_host(i, c) - _domain_min(c) + _period(c), _period(c));
    const auto * elem = (*_locator)(Point(x_host(i, 0), x_host(i, 1), x_host(i, 2)));
    if (!elem)
      elem_host(i) = DEM::ParticleCloud::exited;
    else if (elem->processor_id() == processor_id())
      elem_host(i) = mesh.getContiguousElementID(elem);
    else if (const auto it = ghost_ids.find(elem); it != ghost_ids.end())
      elem_host(i) = it->second;
    else
    {
      target_host(i) = elem->processor_id();
      far_host(i) = elem->id();
    }
  }
  if (_periodic.size())
    ::Kokkos::deep_copy(_cloud.x, x_host);
  ::Kokkos::deep_copy(_cloud.elem, elem_host);
  ::Kokkos::deep_copy(_cloud.target_rank, target_host);
  ::Kokkos::deep_copy(_far_elem, far_host);
}

bool
KokkosParticleCloud::updateNeighborList()
{
  // Both sides of a cross-rank pair must rebuild together
  unsigned int stale = _force_rebuild || _neighbor_list.stale(_cloud, _cloud.size());
  _force_rebuild = false;
  if (n_processors() > 1)
    _communicator.max(stale);
  if (!stale)
    return false;

  // A rebuild replaces the ghosts, so a pending update must land first
  finishGhostUpdate();
  exchangeGhosts();

  // Ghosts lie within the cutoff of our box; particles that walked out of the mesh are clamped
  // into the edge bins
  const Real cutoff = 2 * _r_max + _skin;
  auto box = _rank_boxes[processor_id()];
  box.first -= Point(cutoff, cutoff, cutoff);
  box.second += Point(cutoff, cutoff, cutoff);
  _neighbor_list.build(
      _cloud, _cloud.size(), _cloud.size() + _num_ghosts, box.min(), box.max(), _skin, _broad_phase);
  ++_num_neighbor_list_builds;
  if (_verify_neighbor_list)
    verifyNeighborList();
  rebuildPairStates();
  return true;
}

bool
KokkosParticleCloud::historyNeeded() const
{
  return _friction.mu > 0 || _friction.mu_r > 0 || _wall_friction.mu > 0 || _wall_friction.mu_r > 0;
}

void
KokkosParticleCloud::rebuildPairStates()
{
  if (!contactEnabled() || !historyNeeded())
    return;
  const auto cloud = _cloud;
  const auto n_local = _cloud.size();
  const auto walls = _walls;
  const auto sidesets = _sideset_walls;
  const Real skin = _skin;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  const auto ghost_self = _ghost_self;
  const auto old = _pair_states;

  // One entry per listed pair from the side that advances it (see updatePairStates()) and per
  // wall within the skin; the map is rebuilt from scratch so pairs no longer listed drop out.
  // The count bounds the entries, so the capacity is rarely exhausted; when it is, the map is
  // regrown and refilled, inserting being idempotent
  std::size_t needed = 0;
  ::Kokkos::parallel_reduce(
      "dem_pair_states_count",
      n_local,
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        for (auto k = offsets(i); k < offsets(i + 1); ++k)
          sum += pairs(k) >= n_local || i < pairs(k);
        sum += walls.n + DEM::SidesetWalls::max_contacts;
      },
      needed);
  DEM::PairStateMap map(2 * needed + 128);
  while (true)
  {
    ::Kokkos::parallel_for(
        "dem_pair_states_rebuild", n_local, KOKKOS_LAMBDA(const std::size_t i) {
          const auto insert = [&](const DEM::PairKey & key)
          {
            const auto entry = old.find(key);
            map.insert(key, old.valid_at(entry) ? old.value_at(entry) : DEM::PairState());
          };
          for (auto k = offsets(i); k < offsets(i + 1); ++k)
          {
            const auto j = pairs(k);
            const bool lower = cloud.gid(i) < cloud.gid(j);
            if (j >= n_local ? (ghost_self(j - n_local) && !lower) : i > j)
              continue;
            insert(DEM::PairKey::pair(cloud.gid(i), cloud.gid(j)));
          }
          const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
          for (std::size_t w = 0; w < walls.n; ++w)
            if (walls.overlap(w, xi, cloud.r(i)) > -skin)
              insert(DEM::PairKey::wall(cloud.gid(i), w));
          DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
          const auto n_wc =
              sidesets.contacts(cloud.elem(i), xi, cloud.r(i), cloud.r(i) + skin, wc, false);
          for (unsigned int k = 0; k < n_wc; ++k)
            insert(DEM::PairKey::wall(cloud.gid(i), walls.n + static_cast<std::size_t>(wc[k].id)));
        });
    if (!map.failed_insert())
      break;
    map = DEM::PairStateMap(2 * map.capacity());
  }
  _pair_states = map;
}

void
KokkosParticleCloud::exchangeGhosts()
{
  constexpr auto record_size = DEM::ParticleCloud::record_size;
  constexpr auto exited = DEM::ParticleCloud::exited;
  const auto cloud = _cloud;
  const auto n = _cloud.size();
  const Real cutoff = 2 * _r_max + _skin;
  _num_ghosts = 0;
  _ghost_send_counts.assign(_neighbor_ranks.size(), 0);
  _ghost_recv.clear();
  // No early return for a rank without neighbors: the exchange below is collective, and a rank
  // can be neighborless (no elements, or an isolated partition) while others are not

  // Flag, per target, the local particles whose shifted image is within the cutoff of that
  // rank's box, then compact the flags into the send index with a scan. The box is shifted the
  // other way rather than the particle
  const auto num_targets = _ghost_targets.size();
  if (_slot.extent(0) < n * num_targets)
    _slot = ::Kokkos::View<std::size_t *>("dem_migrate_slot", n * num_targets);
  const auto flags = _slot;
  ::Kokkos::View<Real * [6], ::Kokkos::LayoutRight> boxes("dem_ghost_boxes", num_targets);
  ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight> shifts("dem_ghost_shifts", num_targets);
  auto boxes_host = ::Kokkos::create_mirror_view(boxes);
  auto shifts_host = ::Kokkos::create_mirror_view(shifts);
  for (const auto t : index_range(_ghost_targets))
  {
    const auto & box = _rank_boxes[_neighbor_ranks[_ghost_targets[t].neighbor]];
    for (const auto c : make_range(3))
    {
      shifts_host(t, c) = _ghost_targets[t].shift(c);
      boxes_host(t, c) = box.min()(c) - cutoff - shifts_host(t, c);
      boxes_host(t, 3 + c) = box.max()(c) + cutoff - shifts_host(t, c);
    }
  }
  ::Kokkos::deep_copy(boxes, boxes_host);
  ::Kokkos::deep_copy(shifts, shifts_host);
  ::Kokkos::parallel_for(
      "dem_ghost_flag", n * num_targets, KOKKOS_LAMBDA(const std::size_t k) {
        const auto t = k / n, i = k % n;
        bool inside = cloud.elem(i) != exited;
        for (unsigned int c = 0; c < 3; ++c)
          inside = inside && cloud.x(i, c) >= boxes(t, c) && cloud.x(i, c) <= boxes(t, 3 + c);
        flags(k) = inside;
      });
  std::size_t total = 0;
  ::Kokkos::parallel_scan(
      "dem_ghost_scan",
      n * num_targets,
      KOKKOS_LAMBDA(const std::size_t k, std::size_t & sum, const bool final) {
        const auto flag = flags(k);
        if (final)
          flags(k) = flag ? sum : std::numeric_limits<std::size_t>::max();
        sum += flag;
      },
      total);
  if (_ghost_send_index.extent(0) < total)
  {
    _ghost_send_index = ::Kokkos::View<std::size_t *>("dem_ghost_send_index", total);
    _ghost_send_shift =
        ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight>("dem_ghost_send_shift", total);
  }
  if (_ghost_buffer.extent(0) < total * record_size)
  {
    _ghost_buffer = ::Kokkos::View<Real *>("dem_ghost_buffer", total * record_size);
    _ghost_buffer_host = ::Kokkos::create_mirror_view(_ghost_buffer);
  }
  const auto send_index = _ghost_send_index;
  const auto send_shift = _ghost_send_shift;
  const auto buffer = _ghost_buffer;
  ::Kokkos::parallel_for(
      "dem_ghost_pack", n * num_targets, KOKKOS_LAMBDA(const std::size_t k) {
        const auto slot = flags(k);
        if (slot == std::numeric_limits<std::size_t>::max())
          return;
        const auto t = k / n, i = k % n;
        send_index(slot) = i;
        for (unsigned int c = 0; c < 3; ++c)
          send_shift(slot, c) = shifts(t, c);
        cloud.pack(i, 0, &buffer(slot * record_size));
        // The record holds the image position
        for (unsigned int c = 0; c < 3; ++c)
          buffer(slot * record_size + 1 + c) += shifts(t, c);
      });

  // The send index is ordered by target and the targets by rank, so the per-rank counts are the
  // slot boundaries
  auto buffer_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{},
      ::Kokkos::subview(buffer, std::make_pair(std::size_t(0), total * record_size)));
  auto flags_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{}, ::Kokkos::subview(flags, std::make_pair(std::size_t(0), n * num_targets)));
  std::map<processor_id_type, std::vector<Real>> to_send;
  for (const auto t : index_range(_ghost_targets))
  {
    const auto b = _ghost_targets[t].neighbor;
    auto & data = to_send[_neighbor_ranks[b]];
    for (const auto i : make_range(n))
    {
      const auto slot = flags_host(t * n + i);
      if (slot == std::numeric_limits<std::size_t>::max())
        continue;
      data.insert(data.end(), &buffer_host(slot * record_size), &buffer_host(slot * record_size) + record_size);
      ++_ghost_send_counts[b];
    }
  }

  std::map<processor_id_type, std::vector<Real>> received;
  Parallel::push_parallel_vector_data(
      _communicator,
      to_send,
      [&received](processor_id_type source, const std::vector<Real> & data)
      { received[source] = data; });

  // Append the ghosts behind the local particles, grouped by source rank in rank order so the
  // forwarding of each substep lands in fixed slots
  for (const auto & [source, data] : received)
    _num_ghosts += data.size() / record_size;
  const auto end = n + _num_ghosts;
  if (end > _cloud.capacity())
    reserve(std::max(end, 2 * _cloud.capacity()));
  if (_ghost_self.extent(0) < _num_ghosts)
    _ghost_self = ::Kokkos::View<bool *>("dem_ghost_self", _num_ghosts);
  auto ghost_self_host = ::Kokkos::create_mirror_view(_ghost_self);
  std::size_t begin = n;
  for (const auto & [source, data] : received)
  {
    const auto count = data.size() / record_size;
    _ghost_recv.push_back({source, begin, count});
    DEM::ParticleCloud::HostMirror host(count);
    for (const auto r : make_range(count))
    {
      host.unpack(data.data() + r * record_size, r);
      host.elem(r) = DEM::ParticleCloud::unresolved;
      ghost_self_host(begin - n + r) = source == processor_id();
    }
    _cloud.copyFrom(host, begin);
    begin += count;
  }
  ::Kokkos::deep_copy(_ghost_self, ghost_self_host);
  if (_ghost_recv_buffer.extent(0) < _num_ghosts * ghost_stride)
  {
    _ghost_recv_buffer =
        ::Kokkos::View<Real *>("dem_ghost_recv_buffer", _num_ghosts * ghost_stride);
    _ghost_recv_buffer_host = ::Kokkos::create_mirror_view(_ghost_recv_buffer);
  }
}

void
KokkosParticleCloud::startGhostUpdate()
{
  if (_neighbor_ranks.empty())
    return;
  ++_num_ghost_forwards;

  // Positions, velocities, and spins only; the rest of a ghost's state is fixed between rebuilds
  constexpr std::size_t stride = ghost_stride;
  const auto cloud = _cloud;
  const auto send_index = _ghost_send_index;
  const auto send_shift = _ghost_send_shift;
  const auto buffer = _ghost_buffer;
  std::size_t total = 0;
  for (const auto count : _ghost_send_counts)
    total += count;
  ::Kokkos::parallel_for(
      "dem_ghost_forward_pack", total, KOKKOS_LAMBDA(const std::size_t slot) {
        const auto i = send_index(slot);
        for (unsigned int c = 0; c < 3; ++c)
        {
          buffer(slot * stride + c) = cloud.x(i, c) + send_shift(slot, c);
          buffer(slot * stride + 3 + c) = cloud.v(i, c);
          buffer(slot * stride + 6 + c) = cloud.omega(i, c);
        }
      });
  ::Kokkos::fence();

  // MPI reads the device buffers directly when the host can access them or the MPI library is
  // GPU-aware, and host mirrors otherwise
  const bool direct = _device_buffers_to_mpi || _gpu_aware_mpi;
  Real * send_ptr = _ghost_buffer.data();
  Real * recv_ptr = _ghost_recv_buffer.data();
  if (!direct)
  {
    ::Kokkos::deep_copy(::Kokkos::subview(_ghost_buffer_host, std::make_pair(std::size_t(0), total * stride)),
                        ::Kokkos::subview(_ghost_buffer, std::make_pair(std::size_t(0), total * stride)));
    send_ptr = _ghost_buffer_host.data();
    recv_ptr = _ghost_recv_buffer_host.data();
  }

  _ghost_requests.clear();
  const auto n = _cloud.size();
  for (const auto & block : _ghost_recv)
  {
    _ghost_requests.emplace_back();
    MPI_Irecv(recv_ptr + (block.begin - n) * stride,
              block.count * stride,
              MPI_DOUBLE,
              block.source,
              _ghost_tag.value(),
              _communicator.get(),
              &_ghost_requests.back());
  }
  std::size_t slot = 0;
  for (const auto b : index_range(_neighbor_ranks))
  {
    if (_ghost_send_counts[b])
    {
      _ghost_requests.emplace_back();
      MPI_Isend(send_ptr + slot * stride,
                _ghost_send_counts[b] * stride,
                MPI_DOUBLE,
                _neighbor_ranks[b],
                _ghost_tag.value(),
                _communicator.get(),
                &_ghost_requests.back());
    }
    slot += _ghost_send_counts[b];
  }
}

void
KokkosParticleCloud::finishGhostUpdate()
{
  if (_ghost_requests.empty())
    return;
  MPI_Waitall(_ghost_requests.size(), _ghost_requests.data(), MPI_STATUSES_IGNORE);
  _ghost_requests.clear();

  constexpr std::size_t stride = ghost_stride;
  const auto n = _cloud.size();
  const auto num_ghosts = _num_ghosts;
  if (!(_device_buffers_to_mpi || _gpu_aware_mpi))
    ::Kokkos::deep_copy(
        ::Kokkos::subview(_ghost_recv_buffer, std::make_pair(std::size_t(0), num_ghosts * stride)),
        ::Kokkos::subview(_ghost_recv_buffer_host, std::make_pair(std::size_t(0), num_ghosts * stride)));
  const auto cloud = _cloud;
  const auto recv = _ghost_recv_buffer;
  ::Kokkos::parallel_for(
      "dem_ghost_forward_unpack", num_ghosts, KOKKOS_LAMBDA(const std::size_t g) {
        for (unsigned int c = 0; c < 3; ++c)
        {
          cloud.x(n + g, c) = recv(g * stride + c);
          cloud.v(n + g, c) = recv(g * stride + 3 + c);
          cloud.omega(n + g, c) = recv(g * stride + 6 + c);
        }
      });
}

void
KokkosParticleCloud::verifyNeighborList()
{
  const auto n = _cloud.size();
  const auto n_total = n + _num_ghosts;
  auto x = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto r = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.r);
  auto gid = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.gid);
  auto offsets =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _neighbor_list.offsets());
  auto pairs = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _neighbor_list.pairs());

  // Pairs with at least one local particle; ghost-ghost pairs are listed but never used
  std::set<std::pair<std::size_t, std::size_t>> listed, expected;
  for (const auto i : make_range(n))
    for (auto k = offsets(i); k < offsets(i + 1); ++k)
      listed.emplace(std::min(i, pairs(k)), std::max(i, pairs(k)));
  for (const auto i : make_range(n))
    for (const auto j : make_range(i + 1, n_total))
    {
      const Point d(x(i, 0) - x(j, 0), x(i, 1) - x(j, 1), x(i, 2) - x(j, 2));
      if (d.norm() < r(i) + r(j) + _skin)
        expected.emplace(i, j);
    }

  for (const auto & [i, j] : expected)
    if (!listed.count({i, j}))
      mooseError("Neighbor list is missing the pair of particles ", gid(i), " and ", gid(j));
  for (const auto & [i, j] : listed)
    if (!expected.count({i, j}))
      mooseError("Neighbor list has the spurious pair of particles ", gid(i), " and ", gid(j));
}

std::size_t
KokkosParticleCloud::migrate()
{
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;
  constexpr auto exited = DEM::ParticleCloud::exited;
  // Each particle travels as a flat record so no TIMPI type registration is needed; the libMesh
  // element ID it carries is mapped back to a contiguous ID by the receiver
  constexpr auto record_size = DEM::ParticleCloud::record_size;

  // Departing and exited particles are found on device so the common case of none costs one scalar
  const auto cloud = _cloud;
  std::size_t n_sent = 0, n_exited = 0;
  ::Kokkos::parallel_reduce(
      "dem_migrate_count",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sent, std::size_t & gone) {
        sent += cloud.target_rank(i) != invalid_proc;
        gone += cloud.elem(i) == exited;
      },
      n_sent,
      n_exited);
  auto n_removed_global = n_sent + n_exited;
  _communicator.sum(n_removed_global);
  if (!n_removed_global)
    return 0;
  _num_exited_local += n_exited;

  // Prefix sums of the keep and departing masks: a kept particle's slot in the compacted cloud,
  // or a departing particle's slot in the send buffer; exited particles get neither
  const auto slot = _slot;
  ::Kokkos::parallel_scan(
      "dem_migrate_scan_keep",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & n_kept_before, const bool final) {
        const bool keep = cloud.target_rank(i) == invalid_proc && cloud.elem(i) != exited;
        if (final && keep)
          slot(i) = n_kept_before;
        n_kept_before += keep;
      });
  ::Kokkos::parallel_scan(
      "dem_migrate_scan_send",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & n_sent_before, const bool final) {
        const bool send = cloud.target_rank(i) != invalid_proc;
        if (final && send)
          slot(i) = n_sent_before;
        n_sent_before += send;
      });

  // Compact kept particles into the scratch cloud and pack departing ones into the send buffer
  const auto scratch = _scratch;
  const auto libmesh_id = _libmesh_id;
  const auto far_elem = _far_elem;
  const auto send_buffer = _send_buffer;
  const auto send_rank = _send_rank;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;
  ::Kokkos::parallel_for(
      "dem_migrate_pack", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.target_rank(i) != invalid_proc)
        {
          const auto j = slot(i);
          const auto elem = cloud.elem(i);
          cloud.pack(i,
                     elem == unresolved ? far_elem(i) : libmesh_id(elem),
                     &send_buffer(j * record_size));
          send_rank(j) = cloud.target_rank(i);
        }
        else if (cloud.elem(i) != exited)
          scratch.copyFrom(cloud, i, slot(i));
      });
  // The histories of the departing particles' contacts travel with them, packed while the old
  // cloud and its neighbor list are still in place
  std::map<processor_id_type, std::vector<Real>> states_to_send;
  packPairStates(states_to_send);

  std::swap(_cloud, _scratch);
  _cloud.n = cloud.size() - n_sent - n_exited;
  // Ghosts were not compacted along; the rebuild after migration replaces them
  _num_ghosts = 0;

  // Only the departing particles cross to the host
  auto send_buffer_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{},
      ::Kokkos::subview(send_buffer, std::make_pair(std::size_t(0), n_sent * record_size)));
  auto send_rank_host = ::Kokkos::create_mirror_view_and_copy(
      ::Kokkos::HostSpace{}, ::Kokkos::subview(send_rank, std::make_pair(std::size_t(0), n_sent)));
  std::map<processor_id_type, std::vector<Real>> to_send;
  for (const auto j : make_range(n_sent))
    to_send[send_rank_host(j)].insert(to_send[send_rank_host(j)].end(),
                                      &send_buffer_host(j * record_size),
                                      &send_buffer_host(j * record_size) + record_size);

  std::vector<Real> received;
  Parallel::push_parallel_vector_data(
      _communicator,
      to_send,
      [&received](processor_id_type, const std::vector<Real> & data)
      { received.insert(received.end(), data.begin(), data.end()); });
  std::vector<Real> states_received;
  Parallel::push_parallel_vector_data(
      _communicator,
      states_to_send,
      [&states_received](processor_id_type, const std::vector<Real> & data)
      { states_received.insert(states_received.end(), data.begin(), data.end()); });
  unpackPairStates(states_received);
  const auto n_received = received.size() / record_size;
  if (!n_received)
    return n_sent;

  // Append the received particles behind the kept ones
  const auto begin = _cloud.size();
  const auto end = begin + n_received;
  if (end > _cloud.capacity())
    reserve(std::max(end, 2 * _cloud.capacity()));

  const auto & mesh = kokkosMesh();
  DEM::ParticleCloud::HostMirror host(n_received);
  for (const auto r : make_range(n_received))
  {
    const Real * record = received.data() + r * record_size;
    host.unpack(record, r);
    const auto * elem = _fe_problem.mesh().getMesh().elem_ptr(record[record_size - 1]);
    mooseAssert(elem->processor_id() == processor_id(), "Received a particle we do not own");
    host.elem(r) = mesh.getContiguousElementID(elem);
  }
  _cloud.copyFrom(host, begin);
  _cloud.n = end;

  return n_sent;
}

void
KokkosParticleCloud::verifyPairStates()
{
  const auto n = _cloud.size();
  const auto n_total = n + _num_ghosts;
  auto x = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto r = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.r);
  auto gid = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.gid);
  auto elem = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.elem);
  auto offsets =
      ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _neighbor_list.offsets());
  auto pairs = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _neighbor_list.pairs());
  auto points = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _walls.points);
  auto normals = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _walls.normals);
  auto states = ::Kokkos::create_mirror(_pair_states);
  ::Kokkos::deep_copy(states, _pair_states);
  libmesh_ignore(n_total);

  // The listed pairs and wall contacts of the local particles, and those of them overlapping
  const auto key_of = [](const DEM::PairKey & key) { return std::make_pair(key.lower, key.upper); };
  std::set<std::pair<int64_t, int64_t>> listed, overlapping;
  for (const auto i : make_range(n))
  {
    if (elem(i) == DEM::ParticleCloud::exited)
      continue;
    const Point xi(x(i, 0), x(i, 1), x(i, 2));
    for (auto k = offsets(i); k < offsets(i + 1); ++k)
    {
      const auto j = pairs(k);
      if (elem(j) == DEM::ParticleCloud::exited)
        continue;
      const auto key = key_of(DEM::PairKey::pair(gid(i), gid(j)));
      listed.insert(key);
      if ((xi - Point(x(j, 0), x(j, 1), x(j, 2))).norm() < r(i) + r(j))
        overlapping.insert(key);
    }
    for (const auto w : make_range(_walls.n))
    {
      const Point p(points(w, 0), points(w, 1), points(w, 2)),
          normal(normals(w, 0), normals(w, 1), normals(w, 2));
      const Real overlap = r(i) - (xi - p) * normal;
      const auto key = key_of(DEM::PairKey::wall(gid(i), w));
      if (overlap > -_skin)
        listed.insert(key);
      if (overlap > 0)
        overlapping.insert(key);
    }
  }

  // Sideset wall contacts are found on device: per particle, the boundaries within the skin
  // (unreduced, the faces given a history) and those touching (reduced), as (particle, face
  // key, touching) records
  if (_sideset_walls.n)
  {
    const auto cloud = _cloud;
    const auto sidesets = _sideset_walls;
    const Real skin = _skin;
    const auto num_walls = _walls.n;
    constexpr auto max_contacts = DEM::SidesetWalls::max_contacts;
    ::Kokkos::View<int64_t * [3]> records("dem_verify_wall_records", n * 2 * max_contacts);
    ::Kokkos::View<unsigned int *> counts("dem_verify_wall_counts", n);
    ::Kokkos::parallel_for(
        "dem_verify_walls", n, KOKKOS_LAMBDA(const std::size_t i) {
          DEM::WallContact wc[max_contacts];
          const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
          unsigned int count = 0;
          for (const bool touching : {false, true})
          {
            const auto found = sidesets.contacts(cloud.elem(i),
                                                 xi,
                                                 cloud.r(i),
                                                 touching ? cloud.r(i) : cloud.r(i) + skin,
                                                 wc,
                                                 touching);
            for (unsigned int k = 0; k < found; ++k, ++count)
            {
              records(i * 2 * max_contacts + count, 0) = cloud.gid(i);
              records(i * 2 * max_contacts + count, 1) =
                  -1 - static_cast<int64_t>(num_walls) - wc[k].id;
              records(i * 2 * max_contacts + count, 2) = touching;
            }
          }
          counts(i) = count;
        });
    auto records_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, records);
    auto counts_h = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, counts);
    for (const auto i : make_range(n))
      for (const auto k : make_range(counts_h(i)))
      {
        const auto key = std::make_pair(records_h(i * 2 * max_contacts + k, 0),
                                        records_h(i * 2 * max_contacts + k, 1));
        listed.insert(key);
        if (records_h(i * 2 * max_contacts + k, 2))
          overlapping.insert(key);
      }
  }

  std::set<std::pair<int64_t, int64_t>> held;
  for (const auto e : make_range(states.capacity()))
  {
    if (!states.valid_at(e))
      continue;
    const auto key = key_of(states.key_at(e));
    held.insert(key);
    if (states.value_at(e).stretched() && !overlapping.count(key))
      mooseError("The contact history of particle ",
                 key.first,
                 " and ",
                 key.second >= 0 ? "particle " : "wall ",
                 key.second >= 0 ? key.second : -1 - key.second,
                 " has a stretched tangential spring but they do not overlap");
  }
  for (const auto & key : listed)
    if (!held.count(key))
      mooseError("The listed contact of particle ",
                 key.first,
                 " and ",
                 key.second >= 0 ? "particle " : "wall ",
                 key.second >= 0 ? key.second : -1 - key.second,
                 " has no history");
}

void
KokkosParticleCloud::packPairStates(std::map<processor_id_type, std::vector<Real>> & to_send)
{
  constexpr auto invalid_proc = libMesh::DofObject::invalid_processor_id;
  if (!contactEnabled() || !historyNeeded())
    return;
  const auto cloud = _cloud;
  const auto n_local = _cloud.size();
  const auto num_walls = _walls.n;
  const auto sidesets = _sideset_walls;
  const Real skin = _skin;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  const auto states = _pair_states;
  const auto sideset_key = KOKKOS_LAMBDA(const std::size_t i, const DEM::WallContact & wc)
  {
    return DEM::PairKey::wall(cloud.gid(i), num_walls + static_cast<std::size_t>(wc.id));
  };
  const auto sideset_contacts =
      KOKKOS_LAMBDA(const std::size_t i, DEM::WallContact * const out)
  {
    return sidesets.contacts(cloud.elem(i),
                             Moose::Kokkos::Real3(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2)),
                             cloud.r(i),
                             cloud.r(i) + skin,
                             out,
                             false);
  };

  // A departing particle's records are its stretched histories: the key and the spring, as
  // Reals (global IDs are exact below 2^53). Count, scan, fill
  const auto slot = _slot;
  const auto packable = KOKKOS_LAMBDA(const std::size_t i, const DEM::PairKey & key)
  {
    const auto entry = states.find(key);
    return states.valid_at(entry) && states.value_at(entry).stretched();
  };
  std::size_t total = 0;
  ::Kokkos::parallel_scan(
      "dem_pair_states_pack_scan",
      n_local,
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum, const bool final) {
        std::size_t count = 0;
        if (cloud.target_rank(i) != invalid_proc)
        {
          for (auto k = offsets(i); k < offsets(i + 1); ++k)
            count += packable(i, DEM::PairKey::pair(cloud.gid(i), cloud.gid(pairs(k))));
          for (std::size_t w = 0; w < num_walls; ++w)
            count += packable(i, DEM::PairKey::wall(cloud.gid(i), w));
          DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
          const auto n_wc = sideset_contacts(i, wc);
          for (unsigned int k = 0; k < n_wc; ++k)
            count += packable(i, sideset_key(i, wc[k]));
        }
        if (final)
          slot(i) = sum;
        sum += count;
      },
      total);
  if (!total)
    return;
  ::Kokkos::View<Real *> buffer("dem_pair_states_pack", total * pair_state_record_size);
  ::Kokkos::View<processor_id_type *> ranks("dem_pair_states_pack_rank", total);
  ::Kokkos::parallel_for(
      "dem_pair_states_pack", n_local, KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.target_rank(i) == invalid_proc)
          return;
        auto r = slot(i);
        const auto pack = [&](const DEM::PairKey & key)
        {
          if (!packable(i, key))
            return;
          const auto & value = states.value_at(states.find(key));
          Real * record = &buffer(r * pair_state_record_size);
          record[0] = key.lower;
          record[1] = key.upper;
          for (unsigned int c = 0; c < 3; ++c)
          {
            record[2 + c] = value.delta_t[c];
            record[5 + c] = value.delta_r[c];
          }
          ranks(r++) = cloud.target_rank(i);
        };
        for (auto k = offsets(i); k < offsets(i + 1); ++k)
          pack(DEM::PairKey::pair(cloud.gid(i), cloud.gid(pairs(k))));
        for (std::size_t w = 0; w < num_walls; ++w)
          pack(DEM::PairKey::wall(cloud.gid(i), w));
        DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
        const auto n_wc = sideset_contacts(i, wc);
        for (unsigned int k = 0; k < n_wc; ++k)
          pack(sideset_key(i, wc[k]));
      });
  auto buffer_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, buffer);
  auto ranks_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, ranks);
  for (const auto r : make_range(total))
    to_send[ranks_host(r)].insert(to_send[ranks_host(r)].end(),
                                  &buffer_host(r * pair_state_record_size),
                                  &buffer_host(r * pair_state_record_size) +
                                      pair_state_record_size);
}

void
KokkosParticleCloud::unpackPairStates(const std::vector<Real> & received)
{
  const auto n = received.size() / pair_state_record_size;
  if (!n)
    return;
  ::Kokkos::View<Real *, ::Kokkos::HostSpace, ::Kokkos::MemoryUnmanaged> host(
      const_cast<Real *>(received.data()), received.size());
  const auto buffer = ::Kokkos::create_mirror_view_and_copy(::Kokkos::DefaultExecutionSpace::memory_space{}, host);
  // Room for the arrivals; a key this rank already holds (its copy of a pair across the
  // partition) is identical and left alone
  if (_pair_states.capacity() < 2 * (_pair_states.size() + n))
    _pair_states.rehash(2 * (_pair_states.size() + n));
  const auto states = _pair_states;
  ::Kokkos::parallel_for(
      "dem_pair_states_unpack", n, KOKKOS_LAMBDA(const std::size_t r) {
        const Real * record = &buffer(r * pair_state_record_size);
        DEM::PairState value;
        for (unsigned int c = 0; c < 3; ++c)
        {
          value.delta_t[c] = record[2 + c];
          value.delta_r[c] = record[5 + c];
        }
        states.insert(DEM::PairKey{static_cast<int64_t>(record[0]), static_cast<int64_t>(record[1])},
                      value);
      });
  if (_pair_states.failed_insert())
    mooseError("Could not insert the received contact histories");
}

void
KokkosParticleCloud::reserve(const std::size_t capacity)
{
  DEM::ParticleCloud cloud;
  cloud.allocate(capacity);
  cloud.copyFrom(_cloud, _cloud.size());
  cloud.n = _cloud.size();
  _cloud = cloud;
  allocateScratch(capacity);
}

void
KokkosParticleCloud::allocateScratch(const std::size_t capacity)
{
  _scratch.allocate(capacity);
  _slot = ::Kokkos::View<std::size_t *>("dem_migrate_slot", capacity);
  // Worst case is every particle departing
  _send_buffer =
      ::Kokkos::View<Real *>("dem_migrate_send", capacity * DEM::ParticleCloud::record_size);
  _send_rank = ::Kokkos::View<processor_id_type *>("dem_migrate_send_rank", capacity);
  _far_elem = ::Kokkos::View<dof_id_type *>("dem_migrate_far_elem", capacity);
}

void
KokkosParticleCloud::count()
{
  const auto cloud = _cloud;
  constexpr auto unresolved = DEM::ParticleCloud::unresolved;

  std::size_t n_unresolved = 0;
  unsigned int max_hops_taken = 0;
  Real ke = 0, ke_rot = 0;
  ::Kokkos::parallel_reduce(
      "dem_count_unresolved",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, std::size_t & sum) {
        sum += cloud.elem(i) == unresolved;
      },
      n_unresolved);
  ::Kokkos::parallel_reduce(
      "dem_count_hops",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, unsigned int & max) {
        max = cloud.hops(i) > max ? cloud.hops(i) : max;
      },
      ::Kokkos::Max<unsigned int>(max_hops_taken));
  ::Kokkos::parallel_reduce(
      "dem_count_ke",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, Real & sum) {
        sum += 0.5 * cloud.m(i) * libMesh::Kokkos::make_vector_ref(cloud.v, i).norm_sq();
      },
      ke);
  ::Kokkos::parallel_reduce(
      "dem_count_ke_rot",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, Real & sum) {
        sum += 0.5 * cloud.inertia(i) * libMesh::Kokkos::make_vector_ref(cloud.omega, i).norm_sq();
      },
      ke_rot);
  // Spin angular momentum; for a sphere the body and space frames coincide
  for (unsigned int c = 0; c < 3; ++c)
  {
    ::Kokkos::parallel_reduce(
        "dem_count_p",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & sum) { sum += cloud.m(i) * cloud.v(i, c); },
        _linear_momentum_vector[c]);
    ::Kokkos::parallel_reduce(
        "dem_count_L",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & sum) {
          sum += cloud.inertia(i) * cloud.omega(i, c);
        },
        _angular_momentum_vector[c]);
  }

  std::size_t num_contacts = 0;
  Real pe = 0, virial = 0;
  if (_contact_model == ContactModel::HERTZ)
  {
    countContacts(_hertz, num_contacts, pe, virial);
    countWallForces(_hertz);
  }
  else
  {
    countContacts(_linear, num_contacts, pe, virial);
    countWallForces(_linear);
  }

  // Bed extent and center of mass
  for (unsigned int c = 0; c < 3; ++c)
  {
    ::Kokkos::parallel_reduce(
        "dem_count_max_extent",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & max) {
          const Real x = cloud.x(i, c) + cloud.r(i);
          max = x > max ? x : max;
        },
        ::Kokkos::Max<Real>(_max_extent[c]));
    ::Kokkos::parallel_reduce(
        "dem_count_min_extent",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & min) {
          const Real x = cloud.x(i, c) - cloud.r(i);
          min = x < min ? x : min;
        },
        ::Kokkos::Min<Real>(_min_extent[c]));
    ::Kokkos::parallel_reduce(
        "dem_count_com",
        cloud.size(),
        KOKKOS_LAMBDA(const std::size_t i, Real & sum) { sum += cloud.m(i) * cloud.x(i, c); },
        _center_of_mass[c]);
  }
  ::Kokkos::parallel_reduce(
      "dem_count_mass",
      cloud.size(),
      KOKKOS_LAMBDA(const std::size_t i, Real & sum) { sum += cloud.m(i); },
      _total_mass);

  // Local values; finalize() reduces them across ranks every step
  _num_particles = cloud.size();
  _num_exited = _num_exited_local;
  _num_unresolved = n_unresolved;
  _max_hops_taken = max_hops_taken;
  _num_neighbor_pairs = _neighbor_list.numPairs();
  _num_contacts = num_contacts;
  _num_pair_states = _pair_states.size();
  _load_imbalance = cloud.size();
  _kinetic_energy = ke;
  _rotational_kinetic_energy = ke_rot;
  _potential_energy = pe;
  _virial = virial;
  _linear_momentum = magnitude(_linear_momentum_vector);
  _angular_momentum = magnitude(_angular_momentum_vector);
}

template <typename Model>
void
KokkosParticleCloud::countWallForces(const Model & contact)
{
  const auto num = numWalls();
  _wall_force_host.assign(3 * num, 0);
  if (!num || !contact.enabled())
    return;
  if (_wall_force.extent(0) < num)
    _wall_force = ::Kokkos::View<Real * [3], ::Kokkos::LayoutRight>("dem_wall_force", num);
  const auto wall_force = _wall_force;
  ::Kokkos::deep_copy(wall_force, 0);

  // The forces the walls exert on the particles, negated, summed per wall from the current
  // histories (not advanced)
  const auto cloud = _cloud;
  const auto eval = evaluator(contact);
  const auto num_walls = _walls.n;
  constexpr auto exited = DEM::ParticleCloud::exited;
  ::Kokkos::parallel_for(
      "dem_wall_forces", cloud.size(), KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited)
          return;
        DEM::Contact c;
        Real f_n;
        Moose::Kokkos::Real3 f_t, tau_r;
        const auto add = [&](const std::size_t w, const Moose::Kokkos::Real3 & f)
        {
          for (unsigned int k = 0; k < 3; ++k)
            ::Kokkos::atomic_add(&wall_force(w, k), -f(k));
        };
        for (std::size_t w = 0; w < num_walls; ++w)
          if (eval.wall(i, w, c))
          {
            eval.forces(c, DEM::PairKey::wall(cloud.gid(i), w), 1, false, f_n, f_t, tau_r);
            add(w, f_n * c.normal + f_t);
          }
        DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
        const auto n_wc = eval.sidesetContacts(i, cloud.r(i), wc);
        for (unsigned int k = 0; k < n_wc; ++k)
        {
          c = eval.sideset(i, wc[k]);
          eval.forces(c, eval.sidesetKey(i, wc[k]), 1, false, f_n, f_t, tau_r);
          add(num_walls + wc[k].boundary, f_n * c.normal + f_t);
        }
      });
  auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, wall_force);
  for (const auto w : make_range(num))
    for (unsigned int k = 0; k < 3; ++k)
      _wall_force_host[3 * w + k] = host(w, k);
}

void
KokkosParticleCloud::profile(const unsigned int direction,
                             const Real lower,
                             const Real upper,
                             const std::size_t num_bins,
                             std::vector<Real> & sums) const
{
  ::Kokkos::View<Real **> device("dem_profile", num_bins, profile_sums);
  if (_contact_model == ContactModel::HERTZ)
    profile(_hertz, direction, lower, upper, num_bins, device);
  else
    profile(_linear, direction, lower, upper, num_bins, device);
  auto host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, device);
  sums.assign(num_bins * profile_sums, 0);
  for (const auto b : make_range(num_bins))
    for (const auto k : make_range(profile_sums))
      sums[b * profile_sums + k] = host(b, k);
}

template <typename Model>
void
KokkosParticleCloud::profile(const Model & contact,
                             const unsigned int direction,
                             const Real lower,
                             const Real upper,
                             const std::size_t num_bins,
                             ::Kokkos::View<Real **> & sums) const
{
  const auto cloud = _cloud;
  const auto n_local = _cloud.size();
  const auto eval = evaluator(contact);
  const bool enabled = contact.enabled();
  const auto num_walls = _walls.n;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  constexpr auto exited = DEM::ParticleCloud::exited;
  const Real width = (upper - lower) / num_bins;
  // Symmetric tensor components in the order xx, yy, zz, xy, xz, yz
  constexpr unsigned int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};

  ::Kokkos::parallel_for(
      "dem_profile", n_local, KOKKOS_LAMBDA(const std::size_t i) {
        if (cloud.elem(i) == exited)
          return;
        const int bin = static_cast<int>(::Kokkos::floor((cloud.x(i, direction) - lower) / width));
        if (bin < 0 || bin >= static_cast<int>(num_bins))
          return;
        const auto add = [&](const unsigned int k, const Real value)
        { ::Kokkos::atomic_add(&sums(bin, k), value); };
        const Real m = cloud.m(i);
        const Moose::Kokkos::Real3 v(cloud.v(i, 0), cloud.v(i, 1), cloud.v(i, 2));
        add(0, 1);
        add(1, m);
        add(2, 4.0 / 3.0 * 3.14159265358979323846 * cloud.r(i) * cloud.r(i) * cloud.r(i));
        for (unsigned int c = 0; c < 3; ++c)
          add(3 + c, m * v(c));
        for (unsigned int k = 0; k < 6; ++k)
          add(6 + k, m * v(ia[k]) * v(ib[k]));
        if (!enabled)
          return;

        // Contact virial: half of each pair to this particle's bin (the pair is visited from
        // both sides), the whole of a wall contact
        DEM::Contact c;
        Real f_n;
        Moose::Kokkos::Real3 f_t, tau_r;
        const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
        for (auto k = offsets(i); k < offsets(i + 1); ++k)
        {
          const auto j = pairs(k);
          if (cloud.elem(j) == exited || !eval.sphere(i, j, c))
            continue;
          const Real sign = cloud.gid(i) < cloud.gid(j) ? 1 : -1;
          eval.forces(
              c, DEM::PairKey::pair(cloud.gid(i), cloud.gid(j)), sign, false, f_n, f_t, tau_r);
          const Moose::Kokkos::Real3 f = f_n * c.normal + f_t;
          const Moose::Kokkos::Real3 r(
              xi(0) - cloud.x(j, 0), xi(1) - cloud.x(j, 1), xi(2) - cloud.x(j, 2));
          for (unsigned int q = 0; q < 6; ++q)
            add(12 + q, 0.5 * r(ia[q]) * f(ib[q]));
        }
        for (std::size_t w = 0; w < num_walls; ++w)
          if (eval.wall(i, w, c))
          {
            eval.forces(c, DEM::PairKey::wall(cloud.gid(i), w), 1, false, f_n, f_t, tau_r);
            const Moose::Kokkos::Real3 f = f_n * c.normal + f_t;
            // The contact point is r_i behind the center along the normal
            const Moose::Kokkos::Real3 r = cloud.r(i) * c.normal;
            for (unsigned int q = 0; q < 6; ++q)
              add(12 + q, r(ia[q]) * f(ib[q]));
          }
        DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
        const auto n_wc = eval.sidesetContacts(i, cloud.r(i), wc);
        for (unsigned int k = 0; k < n_wc; ++k)
        {
          c = eval.sideset(i, wc[k]);
          eval.forces(c, eval.sidesetKey(i, wc[k]), 1, false, f_n, f_t, tau_r);
          const Moose::Kokkos::Real3 f = f_n * c.normal + f_t;
          const Moose::Kokkos::Real3 r = xi - wc[k].point;
          for (unsigned int q = 0; q < 6; ++q)
            add(12 + q, r(ia[q]) * f(ib[q]));
        }
      });
}

RealVectorValue
KokkosParticleCloud::wallForce(const std::size_t wall) const
{
  mooseAssert(wall < numWalls(), "Wall index out of range");
  return RealVectorValue(
      _wall_force_host[3 * wall], _wall_force_host[3 * wall + 1], _wall_force_host[3 * wall + 2]);
}

std::size_t
KokkosParticleCloud::wallIndex(const BoundaryName & boundary) const
{
  const auto it = std::find(_wall_boundaries.begin(), _wall_boundaries.end(), boundary);
  if (it == _wall_boundaries.end())
    mooseError("Boundary '", boundary, "' is not one of the wall_boundaries of ", name());
  return _walls.n + std::distance(_wall_boundaries.begin(), it);
}

template <typename Model>
void
KokkosParticleCloud::countContacts(const Model & contact,
                                   std::size_t & num_contacts,
                                   Real & pe,
                                   Real & virial) const
{
  // Overlapping pairs, their elastic energy, and their virial r_ij . f_ij, each pair counted once
  // across ranks with the neighbor list's rule; wall contacts add elastic energy only, being
  // external to the particle system
  const auto cloud = _cloud;
  const auto n_local = cloud.size();
  const auto walls = _walls;
  const auto sidesets = _sideset_walls;
  const auto offsets = _neighbor_list.offsets();
  const auto pairs = _neighbor_list.pairs();
  constexpr auto exited = DEM::ParticleCloud::exited;
  if (contact.enabled())
    ::Kokkos::parallel_reduce(
        "dem_count_contacts",
        n_local,
        KOKKOS_LAMBDA(const std::size_t i, std::size_t & contacts, Real & energy, Real & w) {
          if (cloud.elem(i) == exited)
            return;
          const Moose::Kokkos::Real3 xi(cloud.x(i, 0), cloud.x(i, 1), cloud.x(i, 2));
          const Moose::Kokkos::Real3 vi(cloud.v(i, 0), cloud.v(i, 1), cloud.v(i, 2));
          for (auto k = offsets(i); k < offsets(i + 1); ++k)
          {
            const auto j = pairs(k);
            if (cloud.elem(j) == exited || !(j < n_local ? i < j : cloud.gid(i) < cloud.gid(j)))
              continue;
            const Moose::Kokkos::Real3 xj(cloud.x(j, 0), cloud.x(j, 1), cloud.x(j, 2));
            const Real distance = (xi - xj).norm();
            const Real overlap = cloud.r(i) + cloud.r(j) - distance;
            if (overlap <= 0)
              continue;
            const Moose::Kokkos::Real3 vj(cloud.v(j, 0), cloud.v(j, 1), cloud.v(j, 2));
            const Real normal_velocity = (vi - vj).dot_product((1.0 / distance) * (xi - xj));
            const Real r_eff = cloud.r(i) * cloud.r(j) / (cloud.r(i) + cloud.r(j));
            const Real m_eff = cloud.m(i) * cloud.m(j) / (cloud.m(i) + cloud.m(j));
            ++contacts;
            energy += contact.energy(overlap, r_eff);
            w += distance * contact.normalForce(overlap, normal_velocity, r_eff, m_eff);
          }
          for (std::size_t wall = 0; wall < walls.n; ++wall)
          {
            const Real overlap = walls.overlap(wall, xi, cloud.r(i));
            if (overlap > 0)
              energy += contact.energy(overlap, cloud.r(i));
          }
          DEM::WallContact wc[DEM::SidesetWalls::max_contacts];
          const auto n_wc = sidesets.contacts(cloud.elem(i), xi, cloud.r(i), cloud.r(i), wc);
          for (unsigned int k = 0; k < n_wc; ++k)
            if (wc[k].overlap > 0)
              energy += contact.energy(wc[k].overlap, cloud.r(i));
        },
        num_contacts,
        pe,
        virial);
}

Real
KokkosParticleCloud::magnitude(const Real (&vector)[3])
{
  return std::sqrt(Utility::pow<2>(vector[0]) + Utility::pow<2>(vector[1]) +
                   Utility::pow<2>(vector[2]));
}

void
KokkosParticleCloud::verify()
{
  const auto num_local = kokkosMesh().getNumLocalElements();
  const auto & ghost_ids = kokkosMesh().getGhostElemIdMapping();
  auto x_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.x);
  auto gid_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.gid);
  auto elem_host = ::Kokkos::create_mirror_view_and_copy(::Kokkos::HostSpace{}, _cloud.elem);

  for (const auto i : make_range(_cloud.size()))
  {
    const Point p(x_host(i, 0), x_host(i, 1), x_host(i, 2));
    const auto * located = (*_locator)(p);
    const auto elem = elem_host(i);

    if (elem == DEM::ParticleCloud::exited)
      mooseError("Particle ", gid_host(i), " at ", p, " exited but was not removed");
    else if (elem == DEM::ParticleCloud::unresolved)
    {
      // The locator fallback only leaves a particle unresolved beyond the ghost layer
      if (located && (located->processor_id() == processor_id() || ghost_ids.count(located)))
        mooseError(
            "Particle ", gid_host(i), " at ", p, " is unresolved but is in element ", located->id());
    }
    else if (elem >= num_local)
      mooseError("Particle ", gid_host(i), " at ", p, " is still in a ghost element after migration");
    else if (!_cid_to_elem[elem]->contains_point(p))
      mooseError("Particle ",
                 gid_host(i),
                 " at ",
                 p,
                 " was assigned element ",
                 _cid_to_elem[elem]->id(),
                 " which does not contain it");
  }
}

void
KokkosParticleCloud::finalize()
{
  _communicator.sum(_num_particles);
  _communicator.sum(_num_exited);
  _communicator.sum(_num_unresolved);
  if (_num_unresolved && !_allow_unresolved && !_initial_finalize)
    mooseError(_num_unresolved,
               " particle(s) could not be located in any element at the end of the step; set "
               "allow_unresolved = true to continue with them integrating without walls");
  _initial_finalize = false;
  _communicator.sum(_num_migrated);
  _communicator.max(_max_hops_taken);
  _communicator.sum(_num_neighbor_pairs);
  _communicator.max(_num_neighbor_list_builds);
  _communicator.max(_num_ghost_forwards);
  _communicator.sum(_num_contacts);
  _communicator.sum(_num_pair_states);
  // Largest local count over the mean
  _communicator.max(_load_imbalance);
  _load_imbalance = _num_particles ? _load_imbalance * n_processors() / _num_particles : 1;
  _communicator.sum(_kinetic_energy);
  _communicator.sum(_rotational_kinetic_energy);
  _communicator.sum(_potential_energy);
  _communicator.sum(_virial);
  for (auto & p : _linear_momentum_vector)
    _communicator.sum(p);
  for (unsigned int c = 0; c < 3; ++c)
  {
    _communicator.max(_max_extent[c]);
    _communicator.min(_min_extent[c]);
    _communicator.sum(_center_of_mass[c]);
  }
  _communicator.sum(_total_mass);
  for (unsigned int c = 0; c < 3; ++c)
    _center_of_mass[c] = _total_mass > 0 ? _center_of_mass[c] / _total_mass : 0;
  _communicator.sum(_wall_force_host);
  for (auto & L : _angular_momentum_vector)
    _communicator.sum(L);
  _linear_momentum = magnitude(_linear_momentum_vector);
  _angular_momentum = magnitude(_angular_momentum_vector);
  // Virial theorem over the mesh dimension d: P V = (2 KE + sum r_ij . f_ij) / d
  const auto dim = _fe_problem.mesh().dimension();
  _virial_pressure = (2 * _kinetic_energy + _virial) / (dim * _domain_volume);
}
