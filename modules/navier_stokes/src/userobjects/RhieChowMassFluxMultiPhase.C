//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "RhieChowMassFluxMultiPhase.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NS.h"
#include "VectorCompositeFunctor.h"
#include "PIMPLE.h"
#include "SIMPLE.h"
#include "PetscVectorReader.h"
#include "LinearSystem.h"
#include "LinearFVBoundaryCondition.h"
#include "BalancedForceSurfaceTension.h"

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/elem_range.h"
#include "libmesh/petsc_matrix.h"

using namespace libMesh;

registerMooseObject("NavierStokesApp", RhieChowMassFluxMultiPhase);

InputParameters
RhieChowMassFluxMultiPhase::validParams()
{
  auto params = RhieChowFaceFluxProvider::validParams();
  params += NonADFunctorInterface::validParams();

  params.addClassDescription("Computes H/A and 1/A together with face mass fluxes for segregated "
                             "momentum-pressure equations using linear systems.");

  params.addRequiredParam<VariableName>(NS::pressure, "The pressure variable.");
  params.addRequiredParam<VariableName>("u", "The x-component of velocity");
  params.addParam<VariableName>("v", "The y-component of velocity");
  params.addParam<VariableName>("w", "The z-component of velocity");
  params.addRequiredParam<std::string>("p_diffusion_kernel",
                                       "The diffusion kernel acting on the pressure.");

  params.addRequiredParam<MooseFunctorName>(NS::density, "Density functor");
  params.addParam<MooseFunctorName>("alpha", 1.0, "Phase fracion");

  // We disable the execution of this, should only provide functions
  // for the SIMPLE executioner
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_NONE);
  exec_enum = {EXEC_NONE};
  params.suppressParameter<ExecFlagEnum>("execute_on");

  // Pressure projection
  params.addParam<MooseEnum>("pressure_projection_method",
                             MooseEnum("standard consistent", "standard"),
                             "The method to use in the pressure projection for Ainv - "
                             "standard (SIMPLE) or consistent (SIMPLEC)");

  params.addParam<UserObjectName>(
      "plic_surface_tension",
      "PLICSurfaceTension user object for integral surface tension. "
      "When provided, the face flux includes a balanced surface tension correction "
      "computed from PLIC geometry, using the same face stencil as the pressure gradient.");

  params.addParam<Real>("rho_1", 0.0, "Density of phase 1 (for consistent mass-momentum transport).");
  params.addParam<Real>("rho_2", 0.0, "Density of phase 2 (for consistent mass-momentum transport).");
  params.addParam<MooseFunctorName>("vof_alpha", "VOF volume fraction variable for consistent "
                                    "mass flux reconstruction (rhoPhi = (alpha*rho1 + (1-alpha)*rho2)*phi). "
                                    "Required when rho_1 and rho_2 are specified.");
  params.addParam<RealVectorValue>("gravity_vector", RealVectorValue(0, 0, 0),
                                   "Gravity vector for p_rgh buoyancy correction. When non-zero, "
                                   "the pressure variable is treated as p_rgh = p - rho*g*h, and a "
                                   "buoyancy flux Ainv*grad(rho*g*h) is added to HbyA. This makes "
                                   "p_rgh ~ 0 at hydrostatic equilibrium, dramatically improving "
                                   "pressure solve conditioning at high density ratios.");

  params.addParam<bool>("check_executioner", true, "Whether to check the type of the executioner");
  params.addParam<std::string>("property_suffix", "A property suffix to add to the RC created objects.");
  params.addParamNamesToGroup("check_executioner property_suffix", "Advanced");

  return params;
}

RhieChowMassFluxMultiPhase::RhieChowMassFluxMultiPhase(const InputParameters & params)
  : RhieChowFaceFluxProvider(params),
    NonADFunctorInterface(this),
    _moose_mesh(UserObject::_subproblem.mesh()),
    _mesh(_moose_mesh.getMesh()),
    _dim(blocksMaxDimension()),
    _p(dynamic_cast<MooseLinearVariableFVReal *>(
        &UserObject::_subproblem.getVariable(0, getParam<VariableName>(NS::pressure)))),
    _vel(_dim, nullptr),
    _prop_suffix(isParamValid("property_suffix") ? "_" + getParam<std::string>("property_suffix") : ""),
    _HbyA_flux(_moose_mesh, blockIDs(), "HbyA_flux" + _prop_suffix),
    _Ainv(_moose_mesh, blockIDs(), "Ainv" + _prop_suffix),
    _face_mass_flux(
        declareRestartableData<FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>>>(
            "face_flux", _moose_mesh, blockIDs(), "face_values")),
    _rho(getFunctor<Real>(NS::density)),
    _alpha(getFunctor<Real>("alpha")),
    _pressure_projection_method(getParam<MooseEnum>("pressure_projection_method")),
    _plic_st(isParamValid("plic_surface_tension")
                 ? &getUserObject<BalancedForceSurfaceTension>("plic_surface_tension")
                 : nullptr),
    _rho_1(getParam<Real>("rho_1")),
    _rho_2(getParam<Real>("rho_2")),
    _vof_alpha(isParamValid("vof_alpha") ? &getFunctor<Real>("vof_alpha") : nullptr),
    _gravity(getParam<RealVectorValue>("gravity_vector")),
    _has_consistent_mass_flux(false),
    _flux_old_captured(false)
{
  if (!_p)
    paramError(NS::pressure, "the pressure must be a MooseLinearVariableFVReal.");
  checkBlocks(*_p);

  std::vector<std::string> vel_names = {"u", "v", "w"};
  for (const auto i : index_range(_vel))
  {
    _vel[i] = dynamic_cast<MooseLinearVariableFVReal *>(
        &UserObject::_subproblem.getVariable(0, getParam<VariableName>(vel_names[i])));

    if (!_vel[i])
      paramError(vel_names[i], "the velocity must be a MOOSELinearVariableFVReal.");
    checkBlocks(*_vel[i]);
  }

  // Register the elemental/face functors which will be queried in the pressure equation
  for (const auto tid : make_range(libMesh::n_threads()))
  {
    UserObject::_subproblem.addFunctor("Ainv" + _prop_suffix, _Ainv, tid);
    UserObject::_subproblem.addFunctor("HbyA" + _prop_suffix, _HbyA_flux, tid);
  }

  if(getParam<bool>("check_executioner"))
    if (!dynamic_cast<SIMPLE *>(getMooseApp().getExecutioner()) &&
        !dynamic_cast<PIMPLE *>(getMooseApp().getExecutioner()))
      mooseError(this->name(),
                " should only be used with a linear segregated thermal-hydraulics solver!");
}

void
RhieChowMassFluxMultiPhase::linkMomentumPressureSystems(
    const std::vector<LinearSystem *> & momentum_systems,
    LinearSystem & pressure_system,
    const std::vector<unsigned int> & momentum_system_numbers)
{
  _momentum_systems = momentum_systems;
  _momentum_system_numbers = momentum_system_numbers;
  _pressure_system = &pressure_system;
  _global_pressure_system_number = _pressure_system->number();

  _momentum_implicit_systems.clear();
  for (auto & system : _momentum_systems)
  {
    _global_momentum_system_numbers.push_back(system->number());
    _momentum_implicit_systems.push_back(dynamic_cast<LinearImplicitSystem *>(&system->system()));
  }

  setupMeshInformation();
}

void
RhieChowMassFluxMultiPhase::meshChanged()
{
  _HbyA_flux.clear();
  _Ainv.clear();
  _face_mass_flux.clear();
  setupMeshInformation();
}

void
RhieChowMassFluxMultiPhase::initialSetup()
{
  // We fetch the pressure diffusion kernel to ensure that the face flux correction
  // is consistent with the pressure discretization in the Poisson equation.
  std::vector<LinearFVFluxKernel *> flux_kernel;
  auto base_query = _fe_problem.theWarehouse()
                        .query()
                        .template condition<AttribThread>(_tid)
                        .template condition<AttribSysNum>(_p->sys().number())
                        .template condition<AttribSystem>("LinearFVFluxKernel")
                        .template condition<AttribName>(getParam<std::string>("p_diffusion_kernel"))
                        .queryInto(flux_kernel);
  if (flux_kernel.size() != 1)
    paramError(
        "p_diffusion_kernel",
        "The kernel with the given name could not be found or multiple instances were identified.");
  _p_diffusion_kernel = dynamic_cast<LinearFVAnisotropicDiffusion *>(flux_kernel[0]);
  if (!_p_diffusion_kernel)
    paramError("p_diffusion_kernel",
               "The provided diffusion kernel should of type LinearFVAnisotropicDiffusion!");
}

void
RhieChowMassFluxMultiPhase::setupMeshInformation()
{
  // We cache the cell volumes into a petsc vector for corrections here so we can use
  // the optimized petsc operations for the normalization
  _cell_volumes = _pressure_system->currentSolution()->zero_clone();
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
    // We have to check this because the variable might not be defined on the given
    // block
    if (hasBlocks(elem_info->subdomain_id()))
    {
      const auto elem_dof = elem_info->dofIndices()[_global_pressure_system_number][0];
      _cell_volumes->set(elem_dof, elem_info->volume() * elem_info->coordFactor());
    }

  _cell_volumes->close();

  _flow_face_info.clear();
  for (auto & fi : _fe_problem.mesh().faceInfo())
    if (hasBlocks(fi->elemPtr()->subdomain_id()) ||
        (fi->neighborPtr() && hasBlocks(fi->neighborPtr()->subdomain_id())))
      _flow_face_info.push_back(fi);
}

void
RhieChowMassFluxMultiPhase::initialize()
{
  for (const auto & pair : _HbyA_flux)
    _HbyA_flux[pair.first] = 0;

  for (const auto & pair : _Ainv)
    _Ainv[pair.first] = 0;
}

void
RhieChowMassFluxMultiPhase::initFaceMassFlux()
{
  using namespace Moose::FV;

  const auto time_arg = Moose::currentState();
  const bool volumetric_pressure = (_vof_alpha != nullptr);

  // Initialize face fluxes from velocity ICs.
  // In VOF mode, _face_mass_flux stores phi (volumetric) directly.
  // In Euler-Euler mode, it stores rho*alpha*phi (mass flux).
  for (auto & fi : _flow_face_info)
  {
    RealVectorValue face_velocity;

    if (_vel[0]->isInternalFace(*fi))
    {
      const auto & elem_info = *fi->elemInfo();
      const auto & neighbor_info = *fi->neighborInfo();

      if (volumetric_pressure)
      {
        // VOF mode: store phi = u · n directly
        for (const auto dim_i : index_range(_vel))
          interpolate(InterpMethod::Average,
                      face_velocity(dim_i),
                      _vel[dim_i]->getElemValue(elem_info, time_arg),
                      _vel[dim_i]->getElemValue(neighbor_info, time_arg),
                      *fi,
                      true);
      }
      else
      {
        // Euler-Euler mode: store rho*alpha*u · n
        Real elem_rho_alpha = _rho(makeElemArg(fi->elemPtr()), time_arg) * _alpha(makeElemArg(fi->elemPtr()), time_arg);
        Real neighbor_rho_alpha = _rho(makeElemArg(fi->neighborPtr()), time_arg) * _alpha(makeElemArg(fi->neighborPtr()), time_arg);

        for (const auto dim_i : index_range(_vel))
          interpolate(InterpMethod::Average,
                      face_velocity(dim_i),
                      _vel[dim_i]->getElemValue(elem_info, time_arg) * elem_rho_alpha,
                      _vel[dim_i]->getElemValue(neighbor_info, time_arg) * neighbor_rho_alpha,
                      *fi,
                      true);
      }
    }
    else
    {
      const Elem * const boundary_elem =
          hasBlocks(fi->elemPtr()->subdomain_id()) ? fi->elemPtr() : fi->neighborPtr();
      const Moose::FaceArg boundary_face{
          fi, Moose::FV::LimiterType::CentralDifference, true, false, boundary_elem, nullptr};

      if (volumetric_pressure)
      {
        for (const auto dim_i : index_range(_vel))
          face_velocity(dim_i) = raw_value((*_vel[dim_i])(boundary_face, time_arg));
      }
      else
      {
        const Real face_rho_alpha = _rho(boundary_face, time_arg) * _alpha(boundary_face, time_arg);
        for (const auto dim_i : index_range(_vel))
          face_velocity(dim_i) = face_rho_alpha * raw_value((*_vel[dim_i])(boundary_face, time_arg));
      }
    }

    _face_mass_flux[fi->id()] = face_velocity * fi->normal();
  }
}

Real
RhieChowMassFluxMultiPhase::getMassFlux(const FaceInfo & fi) const
{
  if (_vof_alpha)
  {
    // In VOF mode, _face_mass_flux stores phi (volumetric). Reconstruct mass flux.
    const Moose::FaceArg face_arg{&fi,
                                  Moose::FV::LimiterType::CentralDifference,
                                  true, false, fi.elemPtr(), nullptr};
    const Real rho_f = _rho(face_arg, Moose::currentState());
    const Real alpha_f = _alpha(face_arg, Moose::currentState());
    return rho_f * alpha_f * _face_mass_flux.evaluate(&fi);
  }
  return _face_mass_flux.evaluate(&fi);
}

Real
RhieChowMassFluxMultiPhase::getUnweightedMassFlux(const FaceInfo & fi) const
{
  if (_vof_alpha)
  {
    // In VOF mode, _face_mass_flux stores phi. Unweighted mass flux = rho * phi.
    const Moose::FaceArg face_arg{&fi,
                                  Moose::FV::LimiterType::CentralDifference,
                                  true, false, fi.elemPtr(), nullptr};
    const Real rho_f = _rho(face_arg, Moose::currentState());
    return rho_f * _face_mass_flux.evaluate(&fi);
  }
  else
  {
    // Euler-Euler mode: _face_mass_flux = rho*alpha*phi, divide by alpha to get rho*phi
    const Moose::FaceArg face_arg{&fi,
                                  Moose::FV::LimiterType::CentralDifference,
                                  true, false, fi.elemPtr(), nullptr};
    const Real face_alpha = _alpha(face_arg, Moose::currentState());
    return libmesh_map_find(_face_mass_flux, fi.id()) / std::max(face_alpha, 1e-42);
  }
}

Real
RhieChowMassFluxMultiPhase::getVolumetricFaceFlux(const FaceInfo & fi) const
{
  if (_vof_alpha)
  {
    // In VOF mode, _face_mass_flux IS phi (volumetric) directly — the pressure
    // equation outputs phi without rho*alpha multiplication.
    return _face_mass_flux.evaluate(&fi);
  }
  else
  {
    // Euler-Euler mode: _face_mass_flux = rho*alpha*phi, divide to get phi
    const Moose::FaceArg face_arg{&fi,
                                  Moose::FV::LimiterType::CentralDifference,
                                  true, false, fi.elemPtr(), nullptr};
    const Real face_rho_alpha = _rho(face_arg, Moose::currentState()) *
                                _alpha(face_arg, Moose::currentState());
    return libmesh_map_find(_face_mass_flux, fi.id()) / std::max(face_rho_alpha, 1e-42);
  }
}

Real
RhieChowMassFluxMultiPhase::getVolumetricFaceFlux(const Moose::FV::InterpMethod m,
                                        const FaceInfo & fi,
                                        const Moose::StateArg & time,
                                        const THREAD_ID /*tid*/,
                                        bool libmesh_dbg_var(subtract_mesh_velocity)) const
{
  mooseAssert(!subtract_mesh_velocity, "RhieChowMassFluxMultiPhase does not support moving meshes yet!");

  if (m != Moose::FV::InterpMethod::RhieChow)
    mooseError("Interpolation methods other than Rhie-Chow are not supported!");
  if (time.state != Moose::currentState().state)
    mooseError("Older interpolation times are not supported!");

  return getVolumetricFaceFlux(fi);
}

std::unique_ptr<NumericVector<Number>> &
RhieChowMassFluxMultiPhase::getCellVolumes()
{
  return _cell_volumes;
}

void
RhieChowMassFluxMultiPhase::setConsistentFaceAlpha(const FaceInfo & fi, Real alpha_f_mules)
{
  _consistent_face_alpha[fi.id()] = alpha_f_mules;
}

void
RhieChowMassFluxMultiPhase::captureOldFlux()
{
  if (_flux_old_captured)
    return;
  for (const auto * fi : _flow_face_info)
    _face_flux_old[fi->id()] = _face_mass_flux.evaluate(fi);
  _flux_old_captured = true;
}

void
RhieChowMassFluxMultiPhase::computeConsistentMassFlux()
{
  for (const auto * fi : _flow_face_info)
  {
    const Real phi_f = getVolumetricFaceFlux(*fi);

    if (_vof_alpha && (_rho_1 > 0.0 || _rho_2 > 0.0))
    {
      // Prefer MULES-limited face alpha (bounded, from alpha solve).
      // Fall back to current VOF alpha only at initialization (before any MULES solve).
      Real alpha_f;
      const auto it = _consistent_face_alpha.find(fi->id());
      if (it != _consistent_face_alpha.end())
        alpha_f = it->second;
      else
      {
        const Moose::FaceArg face_arg{fi,
                                      Moose::FV::LimiterType::CentralDifference,
                                      true, false, fi->elemPtr(), nullptr};
        alpha_f = (*_vof_alpha)(face_arg, Moose::currentState());
      }

      // rhoPhi = (alpha*rho1 + (1-alpha)*rho2) * phi  [OpenFOAM interFoam formula]
      const Real rho_mix_f = alpha_f * _rho_1 + (1.0 - alpha_f) * _rho_2;
      _consistent_mass_flux[fi->id()] = rho_mix_f * phi_f;
    }
    else
    {
      // Euler-Euler mode: mass flux is already stored as rho*alpha*phi
      _consistent_mass_flux[fi->id()] = _face_mass_flux.evaluate(fi);
    }
  }
  _has_consistent_mass_flux = true;
}

Real
RhieChowMassFluxMultiPhase::getConsistentMassFlux(const FaceInfo & fi) const
{
  // Always use precomputed consistent flux — never fall through to central-differenced
  // alpha reconstruction. computeConsistentMassFlux() must have been called before this.
  // This guarantees rhoPhi uses MULES-bounded face alpha (or the initial alpha at startup),
  // matching OpenFOAM where rhoPhi is only ever computed from alphaPhi.
  const auto it = _consistent_mass_flux.find(fi.id());
  if (it != _consistent_mass_flux.end())
    return it->second;

  // If we get here, the face wasn't in _flow_face_info (boundary face not covered
  // by computeConsistentMassFlux). Fall back to getMassFlux.
  return getMassFlux(fi);
}

void
RhieChowMassFluxMultiPhase::computeFaceMassFlux()
{
  using namespace Moose::FV;

  // Reset the old flux flag so captureOldFlux works on the next timestep
  _flux_old_captured = false;

  const auto time_arg = Moose::currentState();

  // Petsc vector reader to make the repeated reading from the vector faster
  PetscVectorReader p_reader(*_pressure_system->system().current_local_solution);

  // We loop through the faces and compute the face fluxes using the pressure gradient
  // and the momentum matrix/right hand side
  for (auto & fi : _flow_face_info)
  {
    // Making sure the kernel knows which face we are on
    _p_diffusion_kernel->setupFaceData(fi);

    // We are setting this to 1.0 because we don't want to multiply the kernel contributions
    // with the surface area yet. The surface area will be factored in in the advection kernels.
    _p_diffusion_kernel->setCurrentFaceArea(1.0);

    Real p_grad_flux = 0.0;
    if (_p->isInternalFace(*fi))
    {
      const auto & elem_info = *fi->elemInfo();
      const auto & neighbor_info = *fi->neighborInfo();

      // Fetching the dof indices for the pressure variable
      const auto elem_dof = elem_info.dofIndices()[_global_pressure_system_number][0];
      const auto neighbor_dof = neighbor_info.dofIndices()[_global_pressure_system_number][0];

      // Fetching the values of the pressure for the element and the neighbor
      auto p_elem_value = p_reader(elem_dof);
      auto p_neighbor_value = p_reader(neighbor_dof);

      // Compute the elem matrix contributions for the face
      const auto elem_matrix_contribution = _p_diffusion_kernel->computeElemMatrixContribution();
      const auto neighbor_matrix_contribution =
          _p_diffusion_kernel->computeNeighborMatrixContribution();
      const auto elem_rhs_contribution =
          _p_diffusion_kernel->computeElemRightHandSideContribution();

      // Compute the face flux from the matrix and right hand side contributions
      p_grad_flux = (p_neighbor_value * neighbor_matrix_contribution +
                     p_elem_value * elem_matrix_contribution) -
                    elem_rhs_contribution;
    }
    else if (auto * bc_pointer = _p->getBoundaryCondition(*fi->boundaryIDs().begin()))
    {
      mooseAssert(fi->boundaryIDs().size() == 1, "We should only have one boundary on every face.");

      bc_pointer->setupFaceData(
          fi, fi->faceType(std::make_pair(_p->number(), _global_pressure_system_number)));

      const ElemInfo & elem_info =
          hasBlocks(fi->elemPtr()->subdomain_id()) ? *fi->elemInfo() : *fi->neighborInfo();
      auto p_elem_value = _p->getElemValue(elem_info, time_arg);

      const auto matrix_contribution =
          _p_diffusion_kernel->computeBoundaryMatrixContribution(*bc_pointer);
      const auto rhs_contribution =
          _p_diffusion_kernel->computeBoundaryRHSContribution(*bc_pointer);

      // On the boundary, only the element side has a contribution
      p_grad_flux = (p_elem_value * matrix_contribution - rhs_contribution);
    }
    // Compute the new face flux.
    // Note: when PLIC surface tension is active, _HbyA_flux already includes
    // the surface tension contribution (added in computeHbyA), using the same
    // discrete diffusion operator as p_grad_flux. At equilibrium:
    //   _HbyA_flux = st_flux (velocity HbyA is zero, only ST remains)
    //   p_grad_flux = st_flux (pressure adjusts to balance ST)
    //   → face_mass_flux = -st_flux + st_flux = 0
    // In VOF mode, this IS phi (volumetric) directly.
    // In Euler-Euler mode, this is rho*alpha*phi (mass flux).
    _face_mass_flux[fi->id()] = -_HbyA_flux[fi->id()] + p_grad_flux;
  }
}

void
RhieChowMassFluxMultiPhase::computeCellVelocity()
{

  // Populating cell alpha
  const auto cell_alpha = _pressure_system->currentSolution()->zero_clone();
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
    if (hasBlocks(elem_info->subdomain_id()))
    {
      const auto elem_dof = elem_info->dofIndices()[_global_pressure_system_number][0];
      cell_alpha->set(elem_dof, _alpha(makeElemArg(elem_info->elem()), Moose::currentState()));
    }
  cell_alpha->close();

  auto & pressure_gradient = _pressure_system->gradientContainer();

  // u_C = -(H/A)_C - (1/A)_C * alpha_C * grad(p)_C where C is the cell index
  // In VOF mode, alpha is 1.0 everywhere (single velocity field), so skip the multiplication
  const bool volumetric_pressure = (_vof_alpha != nullptr);

  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    auto working_vector = _Ainv_raw[system_i]->clone();
    if (!volumetric_pressure)
      working_vector->pointwise_mult(*working_vector, *cell_alpha);
    working_vector->pointwise_mult(*working_vector, *pressure_gradient[system_i]);
    working_vector->add(*_HbyA_raw[system_i]);
    working_vector->scale(-1.0);
    (*_momentum_implicit_systems[system_i]->solution) = *working_vector;
    _momentum_implicit_systems[system_i]->update();
    _momentum_systems[system_i]->setSolution(
        *_momentum_implicit_systems[system_i]->current_local_solution);
  }
}

void
RhieChowMassFluxMultiPhase::initCouplingField()
{
  // We loop through the faces and populate the coupling fields (face H/A and 1/H)
  // with 0s for now. Pressure corrector solves will always come after the
  // momentum source so we expect these fields to change before the actual solve.
  for (auto & fi : _fe_problem.mesh().faceInfo())
  {
    _Ainv[fi->id()];
    _HbyA_flux[fi->id()];
  }
}

void
RhieChowMassFluxMultiPhase::populateCouplingFunctors(
    const std::vector<std::unique_ptr<NumericVector<Number>>> & raw_hbya,
    const std::vector<std::unique_ptr<NumericVector<Number>>> & raw_Ainv)
{
  // We have the raw H/A and 1/A vectors in a petsc format. This function
  // will create face functors from them
  using namespace Moose::FV;
  const auto time_arg = Moose::currentState();

  // Create the petsc vector readers for faster repeated access
  std::vector<PetscVectorReader> hbya_reader;
  for (const auto dim_i : index_range(raw_hbya))
    hbya_reader.emplace_back(*raw_hbya[dim_i]);

  std::vector<PetscVectorReader> ainv_reader;
  for (const auto dim_i : index_range(raw_Ainv))
    ainv_reader.emplace_back(*raw_Ainv[dim_i]);

  // We loop through the faces and populate the coupling fields (face H/A and 1/H)
  for (auto & fi : _flow_face_info)
  {
    Real face_rho = 0;
    Real face_alpha = 0;
    RealVectorValue face_hbya;

    // We do the lookup in advance
    auto & Ainv = _Ainv[fi->id()];

    // In VOF mode (_vof_alpha set), the pressure equation outputs VOLUMETRIC flux phi
    // directly, matching OpenFOAM interFoam. Since A() from the momentum equation already
    // contains rho (from ddt(rho,U)), 1/A and H/A are already volumetric — no rho*alpha
    // multiplication needed. In Euler-Euler mode, rho*alpha scaling is required.
    const bool volumetric_pressure = (_vof_alpha != nullptr);

    // If it is internal, we just interpolate (using geometric weights) to the face
    if (_vel[0]->isInternalFace(*fi))
    {
      // Get the dof indices for the element and the neighbor
      const auto & elem_info = *fi->elemInfo();
      const auto & neighbor_info = *fi->neighborInfo();
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[0]][0];
      const auto neighbor_dof = neighbor_info.dofIndices()[_global_momentum_system_numbers[0]][0];

      Real elem_scale = 1.0, neighbor_scale = 1.0;
      if (!volumetric_pressure)
      {
        elem_scale = _rho(makeElemArg(fi->elemPtr()), time_arg) *
                     _alpha(makeElemArg(fi->elemPtr()), time_arg);
        neighbor_scale = _rho(makeElemArg(fi->neighborPtr()), time_arg) *
                         _alpha(makeElemArg(fi->neighborPtr()), time_arg);
        interpolate(Moose::FV::InterpMethod::Average, face_rho, elem_scale, neighbor_scale, *fi, true);
        interpolate(Moose::FV::InterpMethod::Average, face_alpha,
                    _alpha(makeElemArg(fi->elemPtr()), time_arg),
                    _alpha(makeElemArg(fi->neighborPtr()), time_arg), *fi, true);
      }

      for (const auto dim_i : index_range(raw_hbya))
      {
        interpolate(Moose::FV::InterpMethod::Average,
                    face_hbya(dim_i),
                    hbya_reader[dim_i](elem_dof),
                    hbya_reader[dim_i](neighbor_dof),
                    *fi,
                    true);
        if (volumetric_pressure)
        {
          // phi mode: Ainv = V/A (no rho*alpha, since A already has rho)
          interpolate(InterpMethod::Average,
                      Ainv(dim_i),
                      ainv_reader[dim_i](elem_dof),
                      ainv_reader[dim_i](neighbor_dof),
                      *fi,
                      true);
        }
        else
        {
          // Mass flux mode: Ainv = rho * alpha^2 * (1/A) * V
          interpolate(InterpMethod::Average,
                      Ainv(dim_i),
                      elem_scale * _alpha(makeElemArg(fi->elemPtr()), time_arg) * ainv_reader[dim_i](elem_dof),
                      neighbor_scale * _alpha(makeElemArg(fi->neighborPtr()), time_arg) * ainv_reader[dim_i](neighbor_dof),
                      *fi,
                      true);
        }
      }
    }
    else
    {
      const ElemInfo & elem_info =
          hasBlocks(fi->elemPtr()->subdomain_id()) ? *fi->elemInfo() : *fi->neighborInfo();
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[0]][0];

      if (_vel[0]->isDirichletBoundaryFace(*fi))
      {
        const Moose::FaceArg boundary_face{
            fi, Moose::FV::LimiterType::CentralDifference, true, false, elem_info.elem(), nullptr};

        if (!volumetric_pressure)
        {
          face_rho = _rho(boundary_face, Moose::currentState());
          face_alpha = _alpha(boundary_face, Moose::currentState());
        }

        for (const auto dim_i : make_range(_dim))
          face_hbya(dim_i) =
              -MetaPhysicL::raw_value((*_vel[dim_i])(boundary_face, Moose::currentState()));
      }
      else
      {
        const auto elem_dof_inner = elem_info.dofIndices()[_global_momentum_system_numbers[0]][0];

        if (!volumetric_pressure)
        {
          face_rho = _rho(makeElemArg(elem_info.elem()), time_arg);
          face_alpha = _alpha(makeElemArg(elem_info.elem()), time_arg);
        }

        for (const auto dim_i : make_range(_dim))
          face_hbya(dim_i) = hbya_reader[dim_i](elem_dof_inner);
      }

      if (volumetric_pressure)
      {
        for (const auto dim_i : index_range(raw_Ainv))
          Ainv(dim_i) = ainv_reader[dim_i](elem_dof);
      }
      else
      {
        const Real elem_rho = _rho(makeElemArg(elem_info.elem()), time_arg);
        const Real elem_alpha = _alpha(makeElemArg(elem_info.elem()), time_arg);
        for (const auto dim_i : index_range(raw_Ainv))
          Ainv(dim_i) = elem_rho * Utility::pow<2>(elem_alpha) * ainv_reader[dim_i](elem_dof);
      }
    }

    // Populate face HbyA flux
    if (volumetric_pressure)
      _HbyA_flux[fi->id()] = face_hbya * fi->normal();  // phi = HbyA · n (volumetric)
    else
      _HbyA_flux[fi->id()] = face_hbya * fi->normal() * face_rho * face_alpha;  // mass flux

    // ddtCorr: temporal flux correction (OpenFOAM interFoam stabilizer).
    // Adds rho_f * rAU_f * (phi_old - flux(U_old)) / dt to the predicted flux.
    // Since rho*rAU ~ dt (uniform), this doesn't amplify the density contrast.
    // For the first timestep or when no old flux is available, this is zero.
    if (volumetric_pressure && !_face_flux_old.empty())
    {
      const auto phi_old_it = _face_flux_old.find(fi->id());
      if (phi_old_it != _face_flux_old.end())
      {
        const Real phi_old = phi_old_it->second;

        // flux(U_old) = interpolated old velocity · n
        Real flux_U_old = 0.0;
        if (_vel[0]->isInternalFace(*fi))
        {
          const auto & elem_info = *fi->elemInfo();
          const auto & neighbor_info = *fi->neighborInfo();
          RealVectorValue face_vel_old;
          for (const auto dim_i : index_range(_vel))
            interpolate(Moose::FV::InterpMethod::Average,
                        face_vel_old(dim_i),
                        _vel[dim_i]->getElemValue(elem_info, time_arg),
                        _vel[dim_i]->getElemValue(neighbor_info, time_arg),
                        *fi,
                        true);
          flux_U_old = face_vel_old * fi->normal();
        }

        // rho_f * rAU_f = rho_f * (dt/rho_f) = dt exactly (uniform across phases).
        // So the correction is simply: phi_old - flux(U_old)
        // (the dt from rho*rAU cancels with the 1/dt from ddtCorr)
        const Real correction = phi_old - flux_U_old;

        // In our sign convention: _HbyA_flux = -phiHbyA, so we subtract
        _HbyA_flux[fi->id()] -= correction;
      }
    }
  }
}

void
RhieChowMassFluxMultiPhase::computeHbyA(const bool with_updated_pressure, bool verbose)
{
  if (verbose)
  {
    _console << "************************************" << std::endl;
    _console << "Computing HbyA" << std::endl;
    _console << "************************************" << std::endl;
  }
  mooseAssert(_momentum_implicit_systems.size() && _momentum_implicit_systems[0],
              "The momentum system shall be linked before calling this function!");

  auto & pressure_gradient = selectPressureGradient(with_updated_pressure);

  _HbyA_raw.clear();
  _Ainv_raw.clear();

    // Populating cell alpha
    const auto cell_alpha = _pressure_system->currentSolution()->zero_clone();
    for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
      if (hasBlocks(elem_info->subdomain_id()))
      {
        const auto elem_dof = elem_info->dofIndices()[_global_pressure_system_number][0];
        cell_alpha->set(elem_dof, _alpha(makeElemArg(elem_info->elem()), Moose::currentState()));
      }
    cell_alpha->close();

  for (auto system_i : index_range(_momentum_systems))
  {
    LinearImplicitSystem * momentum_system = _momentum_implicit_systems[system_i];

    NumericVector<Number> & rhs = *(momentum_system->rhs);
    NumericVector<Number> & current_local_solution = *(momentum_system->current_local_solution);
    NumericVector<Number> & solution = *(momentum_system->solution);
    PetscMatrix<Number> * mmat = dynamic_cast<PetscMatrix<Number> *>(momentum_system->matrix);
    mooseAssert(mmat,
                "The matrices used in the segregated INSFVRhieChow objects need to be convertable "
                "to PetscMatrix!");

    if (verbose)
    {
      _console << "Matrix in rc object" << std::endl;
      mmat->print();
    }

    // First, we extract the diagonal and we will hold on to it for a little while
    _Ainv_raw.push_back(current_local_solution.zero_clone());
    NumericVector<Number> & Ainv = *(_Ainv_raw.back());

    mmat->get_diagonal(Ainv);

    if (verbose)
    {
      _console << "Velocity solution in H(u)" << std::endl;
      solution.print();
    }

    // Time to create H(u) = M_{offdiag} * u - b_{nonpressure}
    _HbyA_raw.push_back(current_local_solution.zero_clone());
    NumericVector<Number> & HbyA = *(_HbyA_raw.back());

    // We start with the matrix product part, we will do
    // M*u - A*u for 2 reasons:
    // 1, We assume A*u petsc operation is faster than setting the matrix diagonal to 0
    // 2, In PISO loops we need to reuse the matrix so we can't just set the diagonals to 0

    // We create a working vector to ease some of the operations, we initialize its values
    // with the current solution values to have something for the A*u term
    auto working_vector = momentum_system->current_local_solution->zero_clone();
    PetscVector<Number> * working_vector_petsc =
        dynamic_cast<PetscVector<Number> *>(working_vector.get());
    mooseAssert(working_vector_petsc,
                "The vectors used in the RhieChowMassFluxMultiPhase objects need to be convertable "
                "to PetscVectors!");

    mmat->vector_mult(HbyA, solution);
    working_vector_petsc->pointwise_mult(Ainv, solution);
    HbyA.add(-1.0, *working_vector_petsc);

    if (verbose)
    {
      _console << " H(u)" << std::endl;
      HbyA.print();
    }

    // We continue by adding the momentum right hand side contributions
    HbyA.add(-1.0, rhs);

    // Unfortunately, the pressure forces are included in the momentum RHS
    // so we have to correct them back
    working_vector_petsc->pointwise_mult(*pressure_gradient[system_i], *_cell_volumes);
    if (!_vof_alpha)
      working_vector_petsc->pointwise_mult(*working_vector_petsc, *cell_alpha);
    HbyA.add(-1.0, *working_vector_petsc);

    if (verbose)
    {
      _console << "total RHS" << std::endl;
      rhs.print();
      _console << "pressure RHS" << std::endl;
      pressure_gradient[system_i]->print();
      _console << " H(u)-rhs-relaxation_source" << std::endl;
      HbyA.print();
    }

    // It is time to create element-wise 1/A-s based on the the diagonal of the momentum matrix
    *working_vector_petsc = 1.0;
    Ainv.pointwise_divide(*working_vector_petsc, Ainv);

    // Create 1/A*(H(u)-RHS)
    HbyA.pointwise_mult(HbyA, Ainv);

    if (verbose)
    {
      _console << " (H(u)-rhs)/A" << std::endl;
      HbyA.print();
    }

    if (_pressure_projection_method == "consistent")
    {

      // Consistent Corrections to SIMPLE
      // 1. Ainv_old = 1/a_p <- Ainv = 1/(a_p + \sum_n a_n)
      // 2. H(u) <- H(u*) + H(u') = H(u*) - (Ainv - Ainv_old) * grad(p) * Vc

      if (verbose)
        _console << "Performing SIMPLEC projection." << std::endl;

      // Lambda function to calculate the sum of diagonal and neighbor coefficients
      auto get_row_sum = [mmat](NumericVector<Number> & sum_vector)
      {
        // Ensure the sum_vector is zeroed out
        sum_vector.zero();

        // Local row size
        const auto local_size = mmat->local_m();

        for (const auto row_i : make_range(local_size))
        {
          // Get all non-zero components of the row of the matrix
          const auto global_index = mmat->row_start() + row_i;
          std::vector<numeric_index_type> indices;
          std::vector<Real> values;
          mmat->get_row(global_index, indices, values);

          // Sum row elements (no absolute values)
          const Real row_sum = std::accumulate(values.cbegin(), values.cend(), 0.0);

          // Add the sum of diagonal and elements to the sum_vector
          sum_vector.add(global_index, row_sum);
        }
        sum_vector.close();
      };

      // Create a temporary vector to store the sum of diagonal and neighbor coefficients
      auto row_sum = current_local_solution.zero_clone();
      get_row_sum(*row_sum);

      // Create vector with new inverse projection matrix
      auto Ainv_full = current_local_solution.zero_clone();
      *working_vector_petsc = 1.0;
      Ainv_full->pointwise_divide(*working_vector_petsc, *row_sum);
      const auto Ainv_full_old = Ainv_full->clone();

      // Correct HbyA
      Ainv_full->add(-1.0, Ainv);
      working_vector_petsc->pointwise_mult(*Ainv_full, *pressure_gradient[system_i]);
      working_vector_petsc->pointwise_mult(*working_vector_petsc, *_cell_volumes);
      HbyA.add(-1.0, *working_vector_petsc);

      // Correct Ainv
      Ainv = *Ainv_full_old;
    }

    Ainv.pointwise_mult(Ainv, *_cell_volumes);
  }

  // We fill the 1/A and H/A functors
  populateCouplingFunctors(_HbyA_raw, _Ainv_raw);

  // Add PLIC surface tension flux to HbyA so the pressure Poisson equation
  // includes the surface tension source term. The face flux is:
  //   phi_f = -HbyA_flux + Ainv*grad(p) - Ainv*sigma*kappa*grad(alpha)
  // For div(phi) = 0:
  //   div(Ainv*grad(p)) = div(HbyA_flux) + div(Ainv*sigma*kappa*grad(alpha))
  // By adding the ST flux to HbyA_flux, the existing pressure equation
  // div(Ainv*grad(p)) = div(HbyA_flux_modified) captures both terms.
  if (_plic_st)
  {
    for (auto & fi : _flow_face_info)
    {
      if (!_p->isInternalFace(*fi))
        continue;

      const Real kappa_f = _plic_st->getFaceCurvature(*fi);
      if (std::abs(kappa_f) < 1e-30)
        continue;

      const Real alpha_elem = _plic_st->getAlpha(fi->elemInfo()->elem()->id());
      const Real alpha_neigh = _plic_st->getAlpha(fi->neighborInfo()->elem()->id());
      const Real sigma = _plic_st->sigma();

      const Real psi_elem = sigma * kappa_f * alpha_elem;
      const Real psi_neigh = sigma * kappa_f * alpha_neigh;

      // Use the same diffusion operator as the pressure gradient
      _p_diffusion_kernel->setupFaceData(fi);
      _p_diffusion_kernel->setCurrentFaceArea(1.0);

      const auto elem_mc = _p_diffusion_kernel->computeElemMatrixContribution();
      const auto neigh_mc = _p_diffusion_kernel->computeNeighborMatrixContribution();

      _HbyA_flux[fi->id()] += psi_neigh * neigh_mc + psi_elem * elem_mc;
    }
  }

  // p_rgh buoyancy correction (OpenFOAM interFoam formulation).
  // When gravity is non-zero, the pressure variable is p_rgh = p - rho*g*h.
  // The buoyancy flux is ghf * Ainv * snGrad(rho), which acts ONLY at the density
  // interface (snGrad(rho) = 0 in uniform-density regions). The body force rho*g
  // remains in the momentum equation and cancels with grad(p) naturally.
  // Using ghf*snGrad(rho) instead of grad(rho*gh) avoids discrete cancellation errors.
  if (_gravity.norm() > 1e-42)
  {
    for (auto & fi : _flow_face_info)
    {
      if (!_p->isInternalFace(*fi))
        continue;

      // rho at cell centers
      const Real rho_elem = _rho(makeElemArg(fi->elemPtr()), Moose::currentState());
      const Real rho_neigh = _rho(makeElemArg(fi->neighborPtr()), Moose::currentState());

      // ghf = gravity · face_centroid (gravity potential at face)
      const Real ghf = _gravity * fi->faceCentroid();

      // Discrete buoyancy flux: ghf * Ainv * snGrad(rho) * Sf
      // Using the pressure diffusion operator: rho_N*neigh_mc + rho_E*elem_mc ≈ Ainv*snGrad(rho)*Sf
      _p_diffusion_kernel->setupFaceData(fi);
      _p_diffusion_kernel->setCurrentFaceArea(1.0);

      const auto elem_mc = _p_diffusion_kernel->computeElemMatrixContribution();
      const auto neigh_mc = _p_diffusion_kernel->computeNeighborMatrixContribution();

      const Real buoyancy_flux = ghf * (rho_neigh * neigh_mc + rho_elem * elem_mc);

      _HbyA_flux[fi->id()] -= buoyancy_flux;
    }
  }

  if (verbose)
  {
    _console << "************************************" << std::endl;
    _console << "DONE Computing HbyA " << std::endl;
    _console << "************************************" << std::endl;
  }
}

std::vector<std::unique_ptr<NumericVector<Number>>> &
RhieChowMassFluxMultiPhase::selectPressureGradient(const bool updated_pressure)
{
  if (updated_pressure)
  {
    _grad_p_current.clear();
    for (const auto & component : _pressure_system->gradientContainer())
      _grad_p_current.push_back(component->clone());
  }

  return _grad_p_current;
}

