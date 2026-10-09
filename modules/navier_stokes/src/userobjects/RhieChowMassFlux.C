//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "RhieChowMassFlux.h"
#include "SubProblem.h"
#include "MooseMesh.h"
#include "NS.h"
#include "VectorCompositeFunctor.h"
#include "PIMPLE.h"
#include "SIMPLE.h"
#include "PetscVectorReader.h"
#include "LinearSystem.h"
#include "LinearFVGradientManager.h"
#include "LinearFVBoundaryCondition.h"
#include "LinearFVAdvectionDiffusionFunctorDirichletBC.h"
#include "LinearFVPressureCorrectionDiffusion.h"
#include "LinearFVMomentumPressure.h"
#include "LinearFVPressureFluxBC.h"
#include "FVReconstructedPressureGradient.h"
#include "FVUtils.h"
#include "MooseUtils.h"
#include "LinearFVElementalKernel.h"
#include "LinearFVAssemblyConsumer.h"
#include "LinearFVKernel.h"
#include "LinearFVBoundaryCondition.h"
#include "LinearFVAnisotropicDiffusion.h"
#include "LinearFVFluxKernel.h"
#include "LinearFVPressureFluxBC.h"
#include "LinearFVPressureSymmetryBC.h"
#include "MooseUtils.h"

#include <algorithm>
#include <cmath>
#include <limits>

// libMesh includes
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/threads.h"

#include <cmath>

registerMooseObject("NavierStokesApp", RhieChowMassFlux);

struct RhieChowMassFlux::MomentumPredictorOperator : public LinearFVAssemblyConsumer
{
  bool assembly_closed = false;
  bool finalized = false;
  Real relaxation_parameter = 1.0;
  bool enforce_diagonal_dominance = false;
  std::unique_ptr<NumericVector<Number>> rhs_raw;
  std::unique_ptr<NumericVector<Number>> pressure_coupling_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> h_diagonal_correction_raw;
  std::unique_ptr<NumericVector<Number>> interior_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_component_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_relax_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_a_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> boundary_dominance_diagonal_raw;
  std::unique_ptr<NumericVector<Number>> offdiag_abs_sum_raw;
  struct BoundaryMatrixContribution
  {
    dof_id_type face_id;
    dof_id_type dof;
    Real contribution;
  };
  std::vector<BoundaryMatrixContribution> boundary_matrix_contributions;
  mutable libMesh::Threads::spin_mutex assembly_mutex;

  bool complete() const;
  void
  addBoundaryMatrixContribution(dof_id_type face_id, dof_id_type dof, Real contribution) override;
};

InputParameters
RhieChowMassFlux::validParams()
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
  params.addParam<std::vector<std::vector<std::string>>>(
      "body_force_kernel_names",
      {},
      "The body force kernel names."
      "this double vector would have size index_x_dim: 'f1x f2x; f1y f2y; f1z f2z'");

  params.addRequiredParam<MooseFunctorName>(NS::density, "Density functor");

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
  params.addParam<MooseEnum>(
      "pressure_diffusion_interpolation",
      MooseEnum("average harmonic", "average"),
      "The face interpolation method for Ainv in the pressure correction diffusion term.");
  params.addParam<bool>(
      "use_cached_momentum_predictor_operator",
      false,
      "Whether to use a cached assembled/relaxed momentum predictor operator as the starting "
      "point for HbyA/Ainv construction instead of rebuilding that base state from the live "
      "system on demand.");
  params.addParam<bool>(
      "split_momentum_predictor_operator",
      false,
      "Whether the assembled momentum predictor operator already excludes explicit pressure-"
      "gradient and body-force terms, so HbyA should use that operator directly instead of "
      "subtracting those terms back out.");
  return params;
}

RhieChowMassFlux::RhieChowMassFlux(const InputParameters & params)
  : RhieChowFaceFluxProvider(params),
    NonADFunctorInterface(this),
    _moose_mesh(UserObject::_subproblem.mesh()),
    _mesh(_moose_mesh.getMesh()),
    _dim(blocksMaxDimension()),
    _p(dynamic_cast<MooseLinearVariableFVReal *>(
        &UserObject::_subproblem.getVariable(0, getParam<VariableName>(NS::pressure)))),
    _vel(_dim, nullptr),
    _HbyA_flux(_moose_mesh, blockIDs(), "HbyA_flux"),
    _phiHbyA_flux(_moose_mesh, blockIDs(), "phiHbyA_flux"),
    _pressure_predictor_flux(_moose_mesh, blockIDs(), "pressure_predictor_flux"),
    _phig_flux(_moose_mesh, blockIDs(), "phig_flux"),
    _Ainv(_moose_mesh, blockIDs(), "Ainv"),
    _pressure_Ainv(_moose_mesh, blockIDs(), "pressure_Ainv"),
    _face_mass_flux(
        declareRestartableData<FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>>>(
            "face_flux", _moose_mesh, blockIDs(), "face_values")),
    _pressure_equation_flux(_moose_mesh, blockIDs(), "pressure_equation_flux"),
    _pressure_boundary_normal_gradient(
        _moose_mesh, blockIDs(), "pressure_boundary_normal_gradient"),
    _boundary_velocity_face_values(_dim),
    _use_cached_momentum_predictor_operator(
        getParam<bool>("use_cached_momentum_predictor_operator")),
    _split_momentum_predictor_operator(getParam<bool>("split_momentum_predictor_operator")),
    _body_force_kernel_names(
        getParam<std::vector<std::vector<std::string>>>("body_force_kernel_names")),
    _rho(getFunctor<Real>(NS::density)),
    _pressure_system(nullptr),
    _pressure_gradient_field(nullptr),
    _base_pressure_gradient_field(nullptr),
    _global_pressure_system_number(0),
    _pressure_projection_method(getParam<MooseEnum>("pressure_projection_method")),
    _pressure_diffusion_interp_method(getParam<MooseEnum>("pressure_diffusion_interpolation") ==
                                              "harmonic"
                                          ? Moose::FV::InterpMethod::HarmonicAverage
                                          : Moose::FV::InterpMethod::Average)
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
  for (const auto tid : make_range(UserObject::_subproblem.numThreads()))
  {
    UserObject::_subproblem.addFunctor("Ainv", _Ainv, tid);
    UserObject::_subproblem.addFunctor("pressure_Ainv", _pressure_Ainv, tid);
    UserObject::_subproblem.addFunctor("HbyA", _HbyA_flux, tid);
    UserObject::_subproblem.addFunctor("phiHbyA", _phiHbyA_flux, tid);
    UserObject::_subproblem.addFunctor("pressure_predictor_flux", _pressure_predictor_flux, tid);
    UserObject::_subproblem.addFunctor("phig", _phig_flux, tid);
    UserObject::_subproblem.addFunctor("pressure_equation_flux", _pressure_equation_flux, tid);
    UserObject::_subproblem.addFunctor(
        "pressure_boundary_normal_gradient", _pressure_boundary_normal_gradient, tid);
  }

  if (!dynamic_cast<SIMPLE *>(getMooseApp().getExecutioner()) &&
      !dynamic_cast<PIMPLE *>(getMooseApp().getExecutioner()))
    mooseError(this->name(),
               " should only be used with a linear segregated thermal-hydraulics solver!");
}

RhieChowMassFlux::~RhieChowMassFlux() = default;

void
RhieChowMassFlux::linkMomentumPressureSystems(
    const std::vector<LinearSystem *> & momentum_systems,
    LinearSystem & pressure_system,
    const std::vector<unsigned int> & momentum_system_numbers)
{
  _momentum_systems = momentum_systems;
  _momentum_system_numbers = momentum_system_numbers;
  _pressure_system = &pressure_system;
  _global_pressure_system_number = _pressure_system->number();

  auto * const pressure_var =
      dynamic_cast<MooseLinearVariableFVReal *>(&_pressure_system->getVariable(0, _p->number()));
  if (!pressure_var)
    mooseError("The pressure variable in system '",
               _pressure_system->name(),
               "' must be a MooseLinearVariableFVReal.");

  // For each momentum system, find the LinearFVMomentumPressure kernels applicable to this RC
  // object, i.e. those whose blocks overlap the blocks this RC object handles. Kernels entirely
  // outside our blocks are irrelevant and ignored; multiple applicable kernels per system are
  // allowed as long as they partition our blocks without overlapping (checked below).
  std::vector<std::vector<const LinearFVMomentumPressure *>> applicable_kernels(
      momentum_systems.size());
  for (const auto system_index : index_range(momentum_systems))
  {
    std::vector<LinearFVElementalKernel *> kernels;
    _fe_problem.theWarehouse()
        .query()
        .template condition<AttribThread>(_tid)
        .template condition<AttribSystem>("LinearFVElementalKernel")
        .template condition<AttribSysNum>(momentum_system_numbers[system_index])
        .queryInto(kernels);

    for (const auto * const kernel : kernels)
    {
      const auto * const candidate = dynamic_cast<const LinearFVMomentumPressure *>(kernel);
      if (candidate && MooseUtils::setsIntersect(candidate->blockIDs(), blockIDs()))
        applicable_kernels[system_index].push_back(candidate);
    }
  }

  for (const auto system_index : index_range(momentum_systems))
  {
    const auto & kernels = applicable_kernels[system_index];
    const auto & system_name = momentum_systems[system_index]->name();

    // Coverage: the applicable kernels must together cover every block this RC object handles.
    std::set<SubdomainID> covered_blocks;
    for (const auto * const kernel : kernels)
      covered_blocks.insert(kernel->blockIDs().begin(), kernel->blockIDs().end());

    std::set<SubdomainID> missing_blocks;
    std::set_difference(blockIDs().begin(),
                        blockIDs().end(),
                        covered_blocks.begin(),
                        covered_blocks.end(),
                        std::inserter(missing_blocks, missing_blocks.end()));

    if (!missing_blocks.empty())
      mooseError("RhieChowMassFlux '",
                 name(),
                 "' found no applicable LinearFVMomentumPressure kernel for momentum system '",
                 system_name,
                 "' on block(s) ",
                 Moose::stringify(missing_blocks),
                 ".");

    // Component and pairwise block-overlap checks among the applicable kernels for this system.
    for (const auto i : index_range(kernels))
    {
      if (kernels[i]->momentumComponent() != system_index)
        mooseError("LinearFVMomentumPressure kernel '",
                   kernels[i]->name(),
                   "' declares momentum component ",
                   kernels[i]->momentumComponent(),
                   " but is linked to momentum system '",
                   system_name,
                   "' (component ",
                   system_index,
                   ") through RhieChowMassFlux '",
                   name(),
                   "'.");

      for (const auto j : make_range(i + 1, kernels.size()))
      {
        std::set<SubdomainID> overlap;
        for (const auto block : kernels[i]->blockIDs())
          if (hasBlocks(block) && kernels[j]->blockIDs().count(block))
            overlap.insert(block);
        if (!overlap.empty())
          mooseError("Momentum system '",
                     system_name,
                     "' has more than one LinearFVMomentumPressure kernel contributing a "
                     "pressure source on the same block(s) for RhieChowMassFlux '",
                     name(),
                     "': '",
                     kernels[i]->name(),
                     "' and '",
                     kernels[j]->name(),
                     "' overlap on block(s) ",
                     Moose::stringify(overlap),
                     ".");
      }
    }
  }

  // Flatten every system's applicable kernels into one list spanning the whole RC object. Even
  // kernels on disjoint blocks/components must share one coupling-gradient field, because H/A
  // has a single field source for this RC object.
  std::vector<const LinearFVMomentumPressure *> momentum_pressure_kernels;
  for (auto & kernels : applicable_kernels)
    momentum_pressure_kernels.insert(
        momentum_pressure_kernels.end(), kernels.begin(), kernels.end());

  const auto * const coupling_kernel = momentum_pressure_kernels.front();
  const auto & coupling_reader = coupling_kernel->pressureGradientField();

  if (dynamic_cast<const FVReconstructedPressureGradient *>(&coupling_reader.method()) &&
      coupling_kernel->pressureVariableNumber() != _p->number())
    mooseError("FVReconstructedPressureGradient can only be used for the pressure variable "
               "registered on RhieChowMassFlux '",
               name(),
               "'.");

  for (const auto * const kernel : momentum_pressure_kernels)
  {
    const auto & reader = kernel->pressureGradientField();

    if (&reader.system() != &coupling_reader.system() ||
        reader.systemNumber() != coupling_reader.systemNumber())
      mooseError("RhieChowMassFlux '",
                 name(),
                 "': momentum pressure kernels '",
                 coupling_kernel->name(),
                 "' (block(s) ",
                 Moose::stringify(coupling_kernel->blockIDs()),
                 ") and '",
                 kernel->name(),
                 "' (block(s) ",
                 Moose::stringify(kernel->blockIDs()),
                 ") must use pressure gradients from the same system.");

    if (reader.variableNumber() != coupling_reader.variableNumber())
      mooseError("RhieChowMassFlux '",
                 name(),
                 "': momentum pressure kernels '",
                 coupling_kernel->name(),
                 "' (pressure variable number ",
                 coupling_reader.variableNumber(),
                 ") and '",
                 kernel->name(),
                 "' (pressure variable number ",
                 reader.variableNumber(),
                 ") must use the same pressure variable.");

    if (&reader.method() != &coupling_reader.method())
      mooseError("RhieChowMassFlux '",
                 name(),
                 "': momentum pressure kernels '",
                 coupling_kernel->name(),
                 "' and '",
                 kernel->name(),
                 "' must use the same pressure gradient method.");

    if (&reader.components() != &coupling_reader.components())
      mooseError("RhieChowMassFlux '",
                 name(),
                 "': momentum pressure kernels '",
                 coupling_kernel->name(),
                 "' and '",
                 kernel->name(),
                 "' must share the same registered pressure gradient field.");
  }

  _pressure_gradient_field = &coupling_reader;

  if (usingReconstructedPressureGradientMethod())
  {
    reconstructedGradientMethod().linkFlowSystem(*this, coupling_reader);
    const auto & base_reader =
        pressure_var->requestCellGradients(reconstructedGradientMethod().baseGradientMethodName());
    _base_pressure_gradient_field = &base_reader;

    mooseAssert(_base_pressure_gradient_field != _pressure_gradient_field,
                "Reconstructed and base pressure gradient readers must be distinct when "
                "FVReconstructedPressureGradient is active.");
  }
  else
    _base_pressure_gradient_field = _pressure_gradient_field;

  _global_momentum_system_numbers.clear();
  _momentum_implicit_systems.clear();
  for (auto & system : _momentum_systems)
  {
    _global_momentum_system_numbers.push_back(system->number());
    _momentum_implicit_systems.push_back(dynamic_cast<LinearImplicitSystem *>(&system->system()));
  }

  if (usingReconstructedPressureGradientMethod())
  {
    checkReconstructedPressureGradientCompatibility();
    if (_p_diffusion_kernel && _p_diffusion_kernel->useNonorthogonalCorrection() &&
        _pressure_projection_method != "consistent")
      paramError("pressure_projection_method",
                 "FVReconstructedPressureGradient with nonorthogonal pressure diffusion requires "
                 "the consistent pressure projection.");
  }

  setupMeshInformation();

  std::set<BoundaryID> velocity_boundary_ids;
  for (const auto dim_i : make_range(_dim))
    for (const auto & [boundary_id, _] : _vel[dim_i]->getBoundaryConditionMap())
      velocity_boundary_ids.insert(boundary_id);

  const auto is_dirichlet =
      [](const MooseLinearVariableFVReal & variable, const BoundaryID boundary_id)
  {
    return dynamic_cast<const LinearFVAdvectionDiffusionFunctorDirichletBC *>(
               variable.getBoundaryCondition(boundary_id)) != nullptr;
  };

  // The legacy boundary HbyA reconstruction uses the x-velocity BC as a proxy for every velocity
  // component, so all components must have the same Dirichlet classification on those boundaries.
  for (const auto boundary_id : velocity_boundary_ids)
  {
    const auto * const pressure_bc = _p->getBoundaryCondition(boundary_id);
    if (dynamic_cast<const LinearFVPressureFluxBC *>(pressure_bc))
      continue;

    const bool velocity_is_dirichlet = is_dirichlet(*_vel[0], boundary_id);
    for (const auto dim_i : make_range(_dim))
      if (is_dirichlet(*_vel[dim_i], boundary_id) != velocity_is_dirichlet)
        mooseError("All velocity components must either have Dirichlet boundary conditions or "
                   "non-Dirichlet boundary conditions on boundary '",
                   _moose_mesh.getBoundaryName(boundary_id),
                   "' when the pressure boundary condition is not a LinearFVPressureFluxBC.");
  }
}

void
RhieChowMassFlux::meshChanged()
{
  _HbyA_flux.clear();
  _phiHbyA_flux.clear();
  _pressure_predictor_flux.clear();
  _phig_flux.clear();
  _Ainv.clear();
  _pressure_Ainv.clear();
  _face_mass_flux.clear();
  _grad_p_current.clear();
  _pressure_equation_flux.clear();
  _pressure_boundary_normal_gradient.clear();
  _pressure_predictor_flux_adjustment.clear();
  _pressure_predictor_base_flux.clear();
  _pressure_equation_flux_valid = false;
  _pressure_boundary_normal_gradient_valid = false;
  _pressure_predictor_face_state_valid = false;
  for (auto & component_face_values : _boundary_velocity_face_values)
    component_face_values.clear();
  _velocity_boundary_state_valid = false;
  clearMomentumPredictorOperatorCache();
  setupMeshInformation();
}

void
RhieChowMassFlux::timestepSetup()
{
  // Nothing to reset when momentum hasn't been linked (e.g. should_solve_momentum = false, used
  // by the restart-freeze pattern).
  if (!_pressure_gradient_field)
    return;

  _grad_p_current.clear();

  if (!usingReconstructedPressureGradientMethod())
    return;

  // Seed the next attempt from the most recently accepted coupling gradient. The linear FV
  // gradient manager restores that state after rejection and checkpoints it for restart/recovery.
  reconstructedGradientMethod().resetForTimeStep(*this);
}

void
RhieChowMassFlux::prepareMomentumPredictor()
{
  mooseAssert(_pressure_gradient_field,
              "The pressure gradient field must be linked before preparing momentum.");

  _grad_p_current.clear();
  for (const auto & component : pressureGradientComponents())
    _grad_p_current.push_back(component->clone());
}

void
RhieChowMassFlux::finalizeMomentumPredictor()
{
  if (usingReconstructedPressureGradientMethod())
    reconstructedGradientMethod().saveLaggedVelocityGradient(*this);
}

void
RhieChowMassFlux::preparePISOCorrector()
{
  if (!usingReconstructedPressureGradientMethod())
    return;

  mooseAssert(!_grad_p_current.empty(),
              "A momentum predictor must be prepared before another PISO corrector.");

  for (const auto component : make_range(dimension()))
    momentumSystem(component).updateFVGradient(velocityVariable(component).requestCellGradients());
  reconstructedGradientMethod().saveLaggedVelocityGradient(*this);
}

void
RhieChowMassFlux::preparePressureRelaxation()
{
  if (!usingReconstructedPressureGradientMethod())
    return;

  // The current pressure solution and corrected face flux are both unrelaxed here. Form the
  // conservative candidate before pressure relaxation replaces that solution.
  reconstructedGradientMethod().computeCandidateFromCorrectedFlux(*this);
  updateCellVelocity(reconstructedGradientMethod().reconstructedCandidate(*this));
}

void
RhieChowMassFlux::finalizePressureCorrector()
{
  if (usingReconstructedPressureGradientMethod())
  {
    // Refresh every ordinary pressure gradient while leaving the old coupling field untouched
    // until the current reconstructed candidate has been relaxed into it.
    _pressure_system->computeGradientsExcept(pressureGradientField());
    reconstructedGradientMethod().finalizeCouplingPressureGradient(
        *this, basePressureGradientComponents());
    _pressure_system->updateFVGradient(pressureGradientField());
  }
  else
  {
    _pressure_system->computeGradients();
    updateCellVelocity(pressureGradientComponents());
  }
}

void
RhieChowMassFlux::initPressureGradient()
{
  _pressure_system->computeGradients();
}

void
RhieChowMassFlux::initialSetup()
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
               "The provided diffusion kernel should be of type LinearFVAnisotropicDiffusion.");

  // We fetch the body forces kernel to ensure that the face flux correction
  // is accurate.

  // Check if components match the dimension.

  if (!_body_force_kernel_names.empty())
  {
    if (_body_force_kernel_names.size() != _dim)
      paramError("body_force_kernel_names",
                 "The dimension of the body force vector does not match the problem dimension.");

    _body_force_kernels.resize(_dim);

    for (const auto dim_i : make_range(_dim))
      for (const auto & force_name : _body_force_kernel_names[dim_i])
      {
        std::vector<LinearFVElementalKernel *> temp_storage;
        auto base_query_force = _fe_problem.theWarehouse()
                                    .query()
                                    .template condition<AttribThread>(_tid)
                                    .template condition<AttribSysNum>(_vel[dim_i]->sys().number())
                                    .template condition<AttribSystem>("LinearFVElementalKernel")
                                    .template condition<AttribName>(force_name)
                                    .queryInto(temp_storage);
        if (temp_storage.size() != 1)
          paramError("body_force_kernel_names",
                     "The kernel with the given name: " + force_name +
                         " could not be found or multiple instances were identified.");
        _body_force_kernels[dim_i].push_back(temp_storage[0]);
      }
  }
}

void
RhieChowMassFlux::setupMeshInformation()
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
    if (isFlowFace(*fi))
      _flow_face_info.push_back(fi);

  for (const auto * fi : _flow_face_info)
  {
    _pressure_boundary_normal_gradient[fi->id()] = 0.0;
    _pressure_predictor_flux_adjustment[fi->id()] = 0.0;
    _pressure_predictor_base_flux[fi->id()] = 0.0;
    _pressure_predictor_flux[fi->id()] = 0.0;
    _pressure_Ainv[fi->id()] = RealVectorValue();
  }
}

template <typename GradientComponent>
void
RhieChowMassFlux::updateCellVelocity(const std::vector<GradientComponent> & pressure_gradient)
{
  // u_C = -(H/A)_C - (1/A)_C * grad(p)_C.
  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    auto working_vector = _Ainv_raw[system_i]->clone();
    working_vector->pointwise_mult(*working_vector, *pressure_gradient[system_i]);
    working_vector->add(*_HbyA_raw[system_i]);
    working_vector->scale(-1.0);
    (*_momentum_implicit_systems[system_i]->solution) = *working_vector;

#ifndef NDEBUG
    {
      // Proof that the assignment above satisfies the momentum-coupling identity
      // u_C + (H/A)_C + Ainv_C * grad(p)_C = 0. The solution was just assigned directly from
      // this formula, so the residual is pure floating-point round-off from the
      // clone/scale/add sequence, not solver error; 1e-10 is a loose absolute bound with
      // headroom over that round-off.
      auto defect = _Ainv_raw[system_i]->clone();
      defect->pointwise_mult(*defect, *pressure_gradient[system_i]);
      defect->add(*_HbyA_raw[system_i]);
      defect->add(*_momentum_implicit_systems[system_i]->solution);
      mooseAssert(MooseUtils::absoluteFuzzyEqual(defect->linfty_norm(), 0.0, 1e-10),
                  "RhieChowMassFlux: cell-velocity update violates the momentum-coupling "
                  "identity u + H/A + Ainv*grad(p) = 0.");
    }
#endif

    _momentum_implicit_systems[system_i]->update();
    _momentum_systems[system_i]->setSolution(
        *_momentum_implicit_systems[system_i]->current_local_solution);
  }
}

void
RhieChowMassFlux::initialize()
{
  for (const auto & pair : _HbyA_flux)
    _HbyA_flux[pair.first] = 0;

  for (const auto & pair : _phiHbyA_flux)
    _phiHbyA_flux[pair.first] = 0;

  for (const auto & pair : _pressure_predictor_flux)
    _pressure_predictor_flux[pair.first] = 0;

  for (const auto & pair : _phig_flux)
    _phig_flux[pair.first] = 0;

  for (const auto & pair : _Ainv)
    _Ainv[pair.first] = 0;

  for (const auto & pair : _pressure_Ainv)
    _pressure_Ainv[pair.first] = 0;

  for (const auto & pair : _pressure_equation_flux)
    _pressure_equation_flux[pair.first] = 0;

  for (const auto & pair : _pressure_boundary_normal_gradient)
    _pressure_boundary_normal_gradient[pair.first] = 0;

  _pressure_predictor_flux_adjustment.clear();
  for (const auto * fi : _flow_face_info)
    _pressure_predictor_flux_adjustment[fi->id()] = 0.0;

  _pressure_equation_flux_valid = false;
  _pressure_boundary_normal_gradient_valid = false;
  _pressure_predictor_face_state_valid = false;
  for (auto & component_face_values : _boundary_velocity_face_values)
    component_face_values.clear();
  _velocity_boundary_state_valid = false;
}

void
RhieChowMassFlux::initFaceMassFlux()
{
  using namespace Moose::FV;

  const auto time_arg = Moose::currentState();

  // We loop through the faces and compute the resulting face fluxes from the
  // initial conditions for velocity
  for (auto & fi : _flow_face_info)
  {
    RealVectorValue density_times_velocity;

    // On internal face we do a regular interpolation with geometric weights
    if (_vel[0]->isInternalFace(*fi))
    {
      const auto & elem_info = *fi->elemInfo();
      const auto & neighbor_info = *fi->neighborInfo();

      Real elem_rho = _rho(makeElemArg(fi->elemPtr()), time_arg);
      Real neighbor_rho = _rho(makeElemArg(fi->neighborPtr()), time_arg);

      for (const auto dim_i : index_range(_vel))
        interpolate(InterpMethod::Average,
                    density_times_velocity(dim_i),
                    _vel[dim_i]->getElemValue(elem_info, time_arg) * elem_rho,
                    _vel[dim_i]->getElemValue(neighbor_info, time_arg) * neighbor_rho,
                    *fi,
                    true);
    }
    // On the boundary, we just take the boundary values
    else
    {
      const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);
      const Moose::FaceArg boundary_face = makeCenteredFaceArg(fi, boundaryElemInfo(fi).elem());

      const Real face_rho = _rho(boundary_face, time_arg);
      for (const auto dim_i : index_range(_vel))
        density_times_velocity(dim_i) =
            boundary_normal_multiplier * face_rho * boundaryVelocityValue(fi, dim_i, time_arg);
    }

    _face_mass_flux[fi->id()] = density_times_velocity * fi->normal();
  }
}

Real
RhieChowMassFlux::getMassFlux(const FaceInfo & fi) const
{
  return _face_mass_flux.evaluate(&fi);
}

Real
RhieChowMassFlux::getVolumetricFaceFlux(const FaceInfo & fi) const
{
  const std::string neighbor_cell_id =
      fi.neighborPtr() ? std::to_string(fi.neighborPtr()->id()) : "none";
  const auto face_flux = _face_mass_flux.find(fi.id());
  if (face_flux == _face_mass_flux.end())
    mooseError("RhieChowMassFlux '",
               name(),
               "' has no corrected mass flux for face ID ",
               fi.id(),
               " adjacent to cell ID ",
               fi.elem().id(),
               " and neighbor cell ID ",
               neighbor_cell_id,
               ".");

  if (!std::isfinite(face_flux->second))
    mooseError("RhieChowMassFlux '",
               name(),
               "' has non-finite corrected mass flux ",
               face_flux->second,
               " on face ID ",
               fi.id(),
               " adjacent to cell ID ",
               fi.elem().id(),
               " and neighbor cell ID ",
               neighbor_cell_id,
               ".");

  const Moose::FaceArg face_arg{&fi,
                                /*limiter_type=*/Moose::FV::LimiterType::CentralDifference,
                                /*elem_is_upwind=*/true,
                                /*correct_skewness=*/false,
                                &fi.elem(),
                                /*state_limiter*/ nullptr};
  const Real face_rho = _rho(face_arg, Moose::currentState());
  if (!std::isfinite(face_rho) || face_rho <= 0.0)
    mooseError("RhieChowMassFlux '",
               name(),
               "' requires finite, positive density for corrected face flux reconstruction, but "
               "density is ",
               face_rho,
               " on face ID ",
               fi.id(),
               " adjacent to cell ID ",
               fi.elem().id(),
               " and neighbor cell ID ",
               neighbor_cell_id,
               ".");

  const Real volumetric_flux = face_flux->second / face_rho;
  if (!std::isfinite(volumetric_flux))
    mooseError("RhieChowMassFlux '",
               name(),
               "' computed non-finite volumetric flux ",
               volumetric_flux,
               " on face ID ",
               fi.id(),
               " adjacent to cell ID ",
               fi.elem().id(),
               " and neighbor cell ID ",
               neighbor_cell_id,
               ".");

  return volumetric_flux;
}

const LinearFVGradientReader &
RhieChowMassFlux::pressureGradientField() const
{
  if (!_pressure_gradient_field)
    mooseError(
        "The pressure gradient field has not been registered for RhieChowMassFlux '", name(), "'.");

  return *_pressure_gradient_field;
}

const LinearFVGradientReader &
RhieChowMassFlux::basePressureGradientField() const
{
  if (!_base_pressure_gradient_field)
    mooseError("The base pressure gradient field has not been registered for RhieChowMassFlux '",
               name(),
               "'.");

  return *_base_pressure_gradient_field;
}

const std::vector<std::unique_ptr<NumericVector<Number>>> &
RhieChowMassFlux::HbyAComponents() const
{
  mooseAssert(!_HbyA_raw.empty(), "HbyA data is not ready.");

  return _HbyA_raw;
}

const std::vector<std::unique_ptr<NumericVector<Number>>> &
RhieChowMassFlux::AinvComponents() const
{
  mooseAssert(!_Ainv_raw.empty(), "Ainv data is not ready.");

  return _Ainv_raw;
}

const std::vector<NumericVector<Number> *> &
RhieChowMassFlux::basePressureGradientComponents() const
{
  return basePressureGradientField().components();
}

unsigned int
RhieChowMassFlux::pressureVariableNumber() const
{
  return _p->number();
}

Real
RhieChowMassFlux::getVolumetricFaceFlux(const Moose::FV::InterpMethod m,
                                        const FaceInfo & fi,
                                        const Moose::StateArg & time,
                                        const THREAD_ID /*tid*/,
                                        bool libmesh_dbg_var(subtract_mesh_velocity)) const
{
  mooseAssert(!subtract_mesh_velocity, "RhieChowMassFlux does not support moving meshes yet!");

  if (m != Moose::FV::InterpMethod::RhieChow)
    mooseError("Interpolation methods other than Rhie-Chow are not supported!");
  if (time.state != Moose::currentState().state)
    mooseError("Older interpolation times are not supported!");

  return getVolumetricFaceFlux(fi);
}

Real
RhieChowMassFlux::maxCourant(const Real dt) const
{
  if (dt <= 0.0)
    return 0.0;

  std::unordered_map<dof_id_type, Real> cell_flux_sum;
  cell_flux_sum.reserve(_flow_face_info.size());

  for (const auto * fi : _flow_face_info)
  {
    if (!fi)
      continue;

    const Real face_measure = faceMeasure(*fi);
    const Real face_volumetric_flux = std::abs(getVolumetricFaceFlux(*fi));
    const Real volumetric_flux = face_volumetric_flux * face_measure;

    if (fi->elemPtr() && hasBlocks(fi->elemPtr()->subdomain_id()))
      cell_flux_sum[fi->elemPtr()->id()] += volumetric_flux;

    if (fi->neighborPtr() && hasBlocks(fi->neighborPtr()->subdomain_id()))
      cell_flux_sum[fi->neighborPtr()->id()] += volumetric_flux;
  }

  Real max_courant = 0.0;
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!hasBlocks(elem_info->subdomain_id()))
      continue;

    const Real volume = elem_info->volume() * elem_info->coordFactor();
    if (volume <= libMesh::TOLERANCE)
      continue;

    const auto it = cell_flux_sum.find(elem_info->elem()->id());
    if (it == cell_flux_sum.end())
      continue;

    const Real courant = 0.5 * dt * it->second / volume;
    max_courant = std::max(max_courant, courant);
  }

  return max_courant;
}

void
RhieChowMassFlux::cachePressureEquationFlux()
{
  // The anisotropic diffusion kernel's nonorthogonal RHS contribution uses
  // cell pressure gradients. When we recache the solved pressure face flux
  // after a pressure solve / relaxation step, refresh those gradients first so
  // the reconstructed pressure-equation flux stays aligned with the current pressure field.
  const_cast<MooseLinearVariableFVReal *>(_p)->computeCellGradients();

  for (auto & fi : _flow_face_info)
    _pressure_equation_flux[fi->id()] = computeDiscretePressureFaceFlux(fi);

  _pressure_equation_flux_valid = true;
}

Real
RhieChowMassFlux::computeDiscretePressureFaceFlux(const FaceInfo * fi) const
{
  using namespace Moose::FV;

  if (!fi)
    return 0.0;

  PetscVectorReader p_reader(*_pressure_system->system().current_local_solution);

  _p_diffusion_kernel->setupFaceData(fi);
  _p_diffusion_kernel->setCurrentFaceArea(pressureDiffusionFaceArea(fi));

  if (_p->isInternalFace(*fi))
  {
    const auto & elem_info = *fi->elemInfo();
    const auto & neighbor_info = *fi->neighborInfo();
    const auto elem_dof = elem_info.dofIndices()[_global_pressure_system_number][0];
    const auto neighbor_dof = neighbor_info.dofIndices()[_global_pressure_system_number][0];
    const auto p_elem_value = p_reader(elem_dof);
    const auto p_neighbor_value = p_reader(neighbor_dof);
    const auto elem_matrix_contribution = _p_diffusion_kernel->computeElemMatrixContribution();
    const auto neighbor_matrix_contribution =
        _p_diffusion_kernel->computeNeighborMatrixContribution();
    const auto elem_rhs_contribution = _p_diffusion_kernel->computeElemRightHandSideContribution();

    return (p_neighbor_value * neighbor_matrix_contribution +
            p_elem_value * elem_matrix_contribution) -
           elem_rhs_contribution;
  }

  if (auto * bc_pointer = pressureBoundaryCondition(fi))
  {
    const ElemInfo & elem_info = boundaryElemInfo(fi);
    const auto elem_dof = elem_info.dofIndices()[_global_pressure_system_number][0];
    const auto p_elem_value = p_reader(elem_dof);
    const auto matrix_contribution =
        _p_diffusion_kernel->computeBoundaryMatrixContribution(*bc_pointer);
    const auto rhs_contribution = _p_diffusion_kernel->computeBoundaryRHSContribution(*bc_pointer);

    return p_elem_value * matrix_contribution - rhs_contribution;
  }

  return 0.0;
}

Real
RhieChowMassFlux::pressurePredictorFlux(const FaceInfo * fi) const
{
  return libmesh_map_find(_HbyA_flux, fi->id());
}

void
RhieChowMassFlux::updatePressurePredictorFaceState()
{
  for (const auto * fi : _flow_face_info)
  {
    const auto face_id = fi->id();
    _phig_flux[face_id] = 0.0;
    _pressure_predictor_base_flux[face_id] = pressurePredictorFlux(fi);
    _phiHbyA_flux[face_id] = _pressure_predictor_base_flux[face_id];
    _pressure_predictor_flux[face_id] = _pressure_predictor_base_flux[face_id];
  }

  _pressure_predictor_face_state_valid = true;
}

void
RhieChowMassFlux::computeFaceMassFlux()
{
  if (!_pressure_predictor_face_state_valid)
    updatePressurePredictorFaceState();

  if (!_pressure_equation_flux_valid)
    cachePressureEquationFlux();

  for (auto & fi : _flow_face_info)
    _face_mass_flux[fi->id()] = -_phiHbyA_flux[fi->id()] + _pressure_equation_flux[fi->id()];
}

void
RhieChowMassFlux::computeCellVelocity()
{
  const auto & pressure_gradient = pressureGradientComponents();

  // We set the dof value in the solution vector the same logic applies:
  // u_C = -(H/A)_C - (1/A)_C*grad(p)_C where C is the cell index
  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    auto working_vector = _Ainv_raw[system_i]->clone();
    working_vector->pointwise_mult(*working_vector, *pressure_gradient[system_i]);
    working_vector->add(*_HbyA_raw[system_i]);
    working_vector->scale(-1.0);
    (*_momentum_implicit_systems[system_i]->solution) = *working_vector;
    _momentum_implicit_systems[system_i]->update();
    _momentum_systems[system_i]->setSolution(
        *_momentum_implicit_systems[system_i]->current_local_solution);
  }

  _velocity_boundary_state_valid = false;
}

Real
RhieChowMassFlux::cellAinvRaw(const unsigned int system_i, const dof_id_type dof) const
{
  mooseAssert(system_i < _Ainv_raw.size() && _Ainv_raw[system_i],
              "Requested Ainv component is not available.");
  return (*_Ainv_raw[system_i])(dof);
}

bool
RhieChowMassFlux::MomentumPredictorOperator::complete() const
{
  return rhs_raw && pressure_coupling_diagonal_raw && h_diagonal_correction_raw &&
         assembly_closed && finalized;
}

void
RhieChowMassFlux::MomentumPredictorOperator::addBoundaryMatrixContribution(
    const dof_id_type face_id, const dof_id_type dof, const Real contribution)
{
  libMesh::Threads::spin_mutex::scoped_lock lock(assembly_mutex);
  boundary_matrix_contributions.push_back({face_id, dof, contribution});
}

const RhieChowMassFlux::MomentumPredictorOperator *
RhieChowMassFlux::cachedMomentumPredictorOperator(const unsigned int system_i) const
{
  if (system_i >= _momentum_predictor_operators.size())
    return nullptr;

  return _momentum_predictor_operators[system_i].get();
}

const NumericVector<Number> *
RhieChowMassFlux::cachedMomentumPredictorDiagonalRaw(const unsigned int system_i) const
{
  const auto * const predictor = cachedMomentumPredictorOperator(system_i);
  return predictor ? predictor->pressure_coupling_diagonal_raw.get() : nullptr;
}

void
RhieChowMassFlux::clearMomentumPredictorOperatorCache()
{
  LinearFVKernel::setAssemblyConsumer(nullptr);
  _momentum_predictor_operators.clear();
  _cached_predictor_operator_valid = false;
  _cached_predictor_operator_finalized = false;
}

void
RhieChowMassFlux::cacheFVSplitMomentumPredictorOperatorAssembly(const unsigned int system_i,
                                                                SparseMatrix<Number> & matrix,
                                                                const NumericVector<Number> & rhs)
{
  if (!_use_cached_momentum_predictor_operator || !_split_momentum_predictor_operator)
    return;

  mooseAssert(system_i < _momentum_predictor_operators.size() &&
                  _momentum_predictor_operators[system_i],
              "A split momentum predictor operator assembly was not started.");

  auto * mmat = dynamic_cast<PetscMatrix<Number> *>(&matrix);
  mooseAssert(mmat, "The split momentum predictor cache requires a PETSc-backed momentum matrix.");

  auto & predictor = *_momentum_predictor_operators[system_i];

  *predictor.rhs_raw = rhs;
  predictor.rhs_raw->close();

  predictor.interior_diagonal_raw->zero();
  predictor.interior_diagonal_raw->close();
  mmat->get_diagonal(*predictor.interior_diagonal_raw);
  predictor.interior_diagonal_raw->close();
  for (const auto & contribution : predictor.boundary_matrix_contributions)
    predictor.interior_diagonal_raw->add(contribution.dof, -contribution.contribution);
  predictor.interior_diagonal_raw->close();

  predictor.offdiag_abs_sum_raw->zero();
  predictor.offdiag_abs_sum_raw->close();

  std::vector<dof_id_type> offdiag_abs_sum_indices;
  std::vector<Real> offdiag_abs_sum_values;
  offdiag_abs_sum_indices.reserve(matrix.row_stop() - matrix.row_start());
  offdiag_abs_sum_values.reserve(matrix.row_stop() - matrix.row_start());

  for (const auto row_dof : make_range(matrix.row_start(), matrix.row_stop()))
  {
    std::vector<numeric_index_type> column_dofs;
    std::vector<Real> coefficients;
    mmat->get_row(row_dof, column_dofs, coefficients);

    Real offdiag_abs_sum = 0.0;
    for (const auto entry_i : index_range(column_dofs))
      if (column_dofs[entry_i] != row_dof)
        offdiag_abs_sum += std::abs(coefficients[entry_i]);

    offdiag_abs_sum_indices.push_back(row_dof);
    offdiag_abs_sum_values.push_back(offdiag_abs_sum);
  }

  predictor.offdiag_abs_sum_raw->insert(offdiag_abs_sum_values, offdiag_abs_sum_indices);
  predictor.offdiag_abs_sum_raw->close();
}

void
RhieChowMassFlux::addMomentumPredictorExplicitForcing(const unsigned int /*system_i*/,
                                                      NumericVector<Number> & /*rhs*/) const
{
}

bool
RhieChowMassFlux::canUseCachedMomentumPredictorOperator()
{
  if (!_use_cached_momentum_predictor_operator)
    return false;

  if (!_split_momentum_predictor_operator)
    return false;

  if (!_cached_predictor_operator_finalized)
    finalizeCachedMomentumPredictorOperators();

  if (!_cached_predictor_operator_valid)
    return false;

  if (_momentum_predictor_operators.size() != _momentum_systems.size())
    return false;

  return std::all_of(_momentum_predictor_operators.begin(),
                     _momentum_predictor_operators.end(),
                     [](const auto & op) { return op && op->complete(); });
}

void
RhieChowMassFlux::updateCachedMomentumPredictorOperatorValidity()
{
  _cached_predictor_operator_valid =
      _momentum_predictor_operators.size() == _momentum_systems.size() &&
      std::all_of(_momentum_predictor_operators.begin(),
                  _momentum_predictor_operators.end(),
                  [](const auto & op) { return op && op->complete(); });
}

void
RhieChowMassFlux::finalizeCachedMomentumPredictorOperators()
{
  if (!_use_cached_momentum_predictor_operator || !_split_momentum_predictor_operator ||
      _cached_predictor_operator_finalized)
    return;

  if (_momentum_predictor_operators.size() != _momentum_systems.size() ||
      std::any_of(_momentum_predictor_operators.begin(),
                  _momentum_predictor_operators.end(),
                  [](const auto & op) { return !op || !op->assembly_closed; }))
    return;

  struct BoundaryFaceCoefficients
  {
    std::vector<Real> coefficient;
    std::vector<dof_id_type> dof;
    std::vector<bool> has_component;
  };

  const auto n_components = _momentum_predictor_operators.size();
  std::unordered_map<dof_id_type, BoundaryFaceCoefficients> boundary_faces;
  for (const auto system_i : index_range(_momentum_predictor_operators))
  {
    auto & predictor = *_momentum_predictor_operators[system_i];
    for (auto * vector : {predictor.boundary_component_diagonal_raw.get(),
                          predictor.boundary_relax_diagonal_raw.get(),
                          predictor.boundary_a_diagonal_raw.get(),
                          predictor.boundary_dominance_diagonal_raw.get()})
    {
      vector->zero();
      vector->close();
    }

    for (const auto & contribution : predictor.boundary_matrix_contributions)
    {
      predictor.boundary_component_diagonal_raw->add(contribution.dof, contribution.contribution);

      auto & face = boundary_faces[contribution.face_id];
      if (face.coefficient.empty())
      {
        face.coefficient.assign(n_components, 0.0);
        face.dof.assign(n_components, DofObject::invalid_id);
        face.has_component.assign(n_components, false);
      }

      face.coefficient[system_i] += contribution.contribution;
      face.dof[system_i] = contribution.dof;
      face.has_component[system_i] = true;
    }
  }

  for (const auto & [face_id, face] : boundary_faces)
  {
    (void)face_id;
    Real component_sum = 0.0;
    Real max_component_magnitude = 0.0;
    Real min_component = std::numeric_limits<Real>::max();
    unsigned int n_present_components = 0;

    for (const auto component_i : index_range(face.coefficient))
      if (face.has_component[component_i])
      {
        const Real coefficient = face.coefficient[component_i];
        component_sum += coefficient;
        max_component_magnitude = std::max(max_component_magnitude, std::abs(coefficient));
        min_component = std::min(min_component, coefficient);
        ++n_present_components;
      }

    if (!n_present_components)
      continue;

    const Real component_average = component_sum / n_present_components;
    for (const auto component_i : index_range(face.coefficient))
      if (face.has_component[component_i])
      {
        const auto dof = face.dof[component_i];
        auto & predictor = *_momentum_predictor_operators[component_i];
        predictor.boundary_dominance_diagonal_raw->add(dof, max_component_magnitude);
        predictor.boundary_relax_diagonal_raw->add(dof, min_component);
        predictor.boundary_a_diagonal_raw->add(dof, component_average);
      }
  }

  for (const auto system_i : index_range(_momentum_predictor_operators))
  {
    auto & predictor = *_momentum_predictor_operators[system_i];

    for (auto * vector : {predictor.boundary_component_diagonal_raw.get(),
                          predictor.boundary_relax_diagonal_raw.get(),
                          predictor.boundary_a_diagonal_raw.get(),
                          predictor.boundary_dominance_diagonal_raw.get()})
      vector->close();

    const auto first_local = predictor.interior_diagonal_raw->first_local_index();
    const auto last_local = predictor.interior_diagonal_raw->last_local_index();
    std::vector<dof_id_type> indices;
    std::vector<Real> pressure_diagonal;
    std::vector<Real> h_diagonal_correction;
    indices.reserve(last_local - first_local);
    pressure_diagonal.reserve(last_local - first_local);
    h_diagonal_correction.reserve(last_local - first_local);

    {
      PetscVectorReader interior_diagonal_reader(*predictor.interior_diagonal_raw);
      PetscVectorReader boundary_component_diagonal_reader(
          *predictor.boundary_component_diagonal_raw);
      PetscVectorReader boundary_relax_diagonal_reader(*predictor.boundary_relax_diagonal_raw);
      PetscVectorReader boundary_a_diagonal_reader(*predictor.boundary_a_diagonal_raw);
      PetscVectorReader boundary_dominance_diagonal_reader(
          *predictor.boundary_dominance_diagonal_raw);
      PetscVectorReader offdiag_abs_sum_reader(*predictor.offdiag_abs_sum_raw);

      for (const auto dof : make_range(first_local, last_local))
      {
        const Real interior_diagonal = interior_diagonal_reader(dof);
        const Real boundary_a_diagonal = boundary_a_diagonal_reader(dof);
        const Real diagonal_for_dominance =
            interior_diagonal + boundary_dominance_diagonal_reader(dof);

        Real predictor_diagonal = diagonal_for_dominance;
        if (predictor.enforce_diagonal_dominance)
          predictor_diagonal =
              std::max(std::abs(diagonal_for_dominance), offdiag_abs_sum_reader(dof));

        predictor_diagonal /= predictor.relaxation_parameter;
        predictor_diagonal -= boundary_relax_diagonal_reader(dof);

        indices.push_back(dof);
        pressure_diagonal.push_back(predictor_diagonal + boundary_a_diagonal);
        h_diagonal_correction.push_back(interior_diagonal - predictor_diagonal -
                                        boundary_component_diagonal_reader(dof) +
                                        boundary_a_diagonal);
      }
    }

    predictor.pressure_coupling_diagonal_raw->zero();
    predictor.pressure_coupling_diagonal_raw->close();
    predictor.pressure_coupling_diagonal_raw->insert(pressure_diagonal, indices);
    predictor.pressure_coupling_diagonal_raw->close();

    predictor.h_diagonal_correction_raw->zero();
    predictor.h_diagonal_correction_raw->close();
    predictor.h_diagonal_correction_raw->insert(h_diagonal_correction, indices);
    predictor.h_diagonal_correction_raw->close();

    predictor.finalized = true;
  }

  _cached_predictor_operator_finalized = true;
  updateCachedMomentumPredictorOperatorValidity();
}

void
RhieChowMassFlux::computePredictorOperatorBase(const unsigned int system_i,
                                               NumericVector<Number> & base_raw,
                                               NumericVector<Number> & diagonal_raw,
                                               const NumericVector<Number> * rhs_override) const
{
  mooseAssert(system_i < _momentum_implicit_systems.size() && _momentum_implicit_systems[system_i],
              "The requested momentum component is not linked to RhieChowMassFlux.");

  auto * momentum_system = _momentum_implicit_systems[system_i];
  auto * mmat = dynamic_cast<PetscMatrix<Number> *>(momentum_system->matrix);
  mooseAssert(mmat,
              "The matrices used in the RhieChowMassFlux objects need to be convertible to "
              "PetscMatrix.");

  const NumericVector<Number> & rhs = rhs_override ? *rhs_override : *(momentum_system->rhs);
  const NumericVector<Number> & solution = *(momentum_system->solution);

  auto diagonal_parallel = solution.zero_clone();
  diagonal_parallel->close();
  mmat->get_diagonal(*diagonal_parallel);
  diagonal_parallel->close();
  diagonal_raw = *diagonal_parallel;
  diagonal_raw.close();

  auto working_vector = solution.zero_clone();
  working_vector->close();
  auto * working_vector_petsc = dynamic_cast<PetscVector<Number> *>(working_vector.get());
  mooseAssert(working_vector_petsc,
              "The vectors used in the RhieChowMassFlux objects need to be convertible to "
              "PetscVectors.");

  auto base_parallel = solution.zero_clone();
  base_parallel->close();
  mmat->vector_mult(*base_parallel, solution);
  base_parallel->close();
  working_vector_petsc->pointwise_mult(*diagonal_parallel, solution);
  working_vector_petsc->close();
  base_parallel->add(-1.0, *working_vector_petsc);
  base_parallel->close();
  base_parallel->add(-1.0, rhs);
  base_parallel->close();
  base_raw = *base_parallel;
  base_raw.close();

  mooseAssert(base_raw.first_local_index() ==
                      momentum_system->current_local_solution->first_local_index() &&
                  base_raw.last_local_index() ==
                      momentum_system->current_local_solution->last_local_index(),
              "Predictor operator cache must use the current-local vector layout.");
}

void
RhieChowMassFlux::beginFVSplitMomentumPredictorOperatorAssembly(const unsigned int system_i)
{
  if (!_use_cached_momentum_predictor_operator)
    return;

  mooseAssert(system_i < _momentum_implicit_systems.size() && _momentum_implicit_systems[system_i],
              "The requested momentum component is not linked to RhieChowMassFlux.");

  if (_momentum_predictor_operators.size() != _momentum_systems.size())
    _momentum_predictor_operators.resize(_momentum_systems.size());

  auto * momentum_system = _momentum_implicit_systems[system_i];
  NumericVector<Number> & current_local_solution = *(momentum_system->current_local_solution);
  auto predictor = std::make_unique<MomentumPredictorOperator>();

  auto initialize_zero_vector =
      [&current_local_solution](std::unique_ptr<NumericVector<Number>> & v)
  {
    v = current_local_solution.zero_clone();
    v->zero();
    v->close();
  };

  for (auto * vector : {&predictor->rhs_raw,
                        &predictor->pressure_coupling_diagonal_raw,
                        &predictor->h_diagonal_correction_raw,
                        &predictor->interior_diagonal_raw,
                        &predictor->boundary_component_diagonal_raw,
                        &predictor->boundary_relax_diagonal_raw,
                        &predictor->boundary_a_diagonal_raw,
                        &predictor->boundary_dominance_diagonal_raw,
                        &predictor->offdiag_abs_sum_raw})
    initialize_zero_vector(*vector);

  LinearFVKernel::setAssemblyConsumer(predictor.get());
  _momentum_predictor_operators[system_i] = std::move(predictor);
  _cached_predictor_operator_valid = false;
}

void
RhieChowMassFlux::completeFVSplitMomentumPredictorOperatorAssembly(
    const unsigned int system_i,
    const Real relaxation_parameter,
    const bool enforce_diagonal_dominance)
{
  LinearFVKernel::setAssemblyConsumer(nullptr);

  if (!_use_cached_momentum_predictor_operator)
    return;

  mooseAssert(system_i < _momentum_predictor_operators.size() &&
                  _momentum_predictor_operators[system_i],
              "A split momentum predictor operator assembly was not started.");
  mooseAssert(relaxation_parameter > 0.0,
              "The momentum equation relaxation parameter must be positive.");

  auto & predictor = *_momentum_predictor_operators[system_i];

  for (auto * vector : {predictor.rhs_raw.get(),
                        predictor.pressure_coupling_diagonal_raw.get(),
                        predictor.h_diagonal_correction_raw.get(),
                        predictor.interior_diagonal_raw.get(),
                        predictor.boundary_component_diagonal_raw.get(),
                        predictor.boundary_relax_diagonal_raw.get(),
                        predictor.boundary_a_diagonal_raw.get(),
                        predictor.boundary_dominance_diagonal_raw.get(),
                        predictor.offdiag_abs_sum_raw.get()})
    vector->close();
  predictor.relaxation_parameter = relaxation_parameter;
  predictor.enforce_diagonal_dominance = enforce_diagonal_dominance;
  predictor.assembly_closed = true;
  predictor.finalized = false;
  _cached_predictor_operator_finalized = false;

  updateCachedMomentumPredictorOperatorValidity();
}

void
RhieChowMassFlux::updateVelocityBoundaryState()
{
  const auto time_arg = Moose::currentState();

  _velocity_boundary_state_valid = false;
  for (auto & component_face_values : _boundary_velocity_face_values)
    component_face_values.clear();

  for (const auto * fi : _flow_face_info)
  {
    if (_vel[0]->isInternalFace(*fi))
      continue;

    RealVectorValue density_times_velocity;
    const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);
    const Moose::FaceArg boundary_face = makeCenteredFaceArg(fi, boundaryElemInfo(fi).elem());
    const Real face_rho = _rho(boundary_face, time_arg);

    for (const auto component : index_range(_vel))
    {
      if (auto * bc_pointer = velocityBoundaryCondition(fi, component))
      {
        _boundary_velocity_face_values[component][fi->id()] = bc_pointer->computeBoundaryValue();
      }
      else
      {
        const ElemInfo & elem_info = boundaryElemInfo(fi);
        _boundary_velocity_face_values[component][fi->id()] =
            _vel[component]->getElemValue(elem_info, time_arg);
      }

      density_times_velocity(component) = boundary_normal_multiplier * face_rho *
                                          _boundary_velocity_face_values[component][fi->id()];
    }

    // Before a pressure correction has been solved, publish the predictor-side
    // boundary phi from the refreshed patch state. But once the pressure
    // equation has produced a corrected face flux, preserve that phi so the
    // boundary face correction stays matched to the cell writeback instead of
    // being overwritten by a later cache refresh.
    if (!_pressure_equation_flux_valid)
      _face_mass_flux[fi->id()] = density_times_velocity * fi->normal();
  }

  _velocity_boundary_state_valid = true;
  postUpdateVelocityBoundaryState();
}

LinearFVBoundaryCondition *
RhieChowMassFlux::velocityBoundaryCondition(const FaceInfo * fi, const unsigned int component) const
{
  mooseAssert(component < _vel.size(), "Velocity component index out of range.");
  return fi ? _vel[component]->getBoundaryCondition(*fi) : nullptr;
}

Real
RhieChowMassFlux::boundaryVelocityValue(const FaceInfo * fi,
                                        const unsigned int component,
                                        const Moose::StateArg & time_arg) const
{
  mooseAssert(fi, "FaceInfo should not be null when evaluating a boundary velocity.");
  mooseAssert(component < _vel.size(), "Velocity component index out of range.");
  mooseAssert(!_vel[component]->isInternalFace(*fi),
              "boundaryVelocityValue should only be called on boundary faces.");

  if (_velocity_boundary_state_valid)
    if (const auto it = _boundary_velocity_face_values[component].find(fi->id());
        it != _boundary_velocity_face_values[component].end())
      return it->second;

  if (auto * bc_pointer = velocityBoundaryCondition(fi, component))
    return bc_pointer->computeBoundaryValue();

  return _vel[component]->getElemValue(boundaryElemInfo(fi), time_arg);
}

Real
RhieChowMassFlux::boundaryMassFluxTarget(const FaceInfo * fi,
                                         const Moose::StateArg & time_arg) const
{
  mooseAssert(fi && !_vel[0]->isInternalFace(*fi),
              "boundaryMassFluxTarget should only be called on boundary faces.");

  const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);
  const Moose::FaceArg boundary_face = makeCenteredFaceArg(fi, boundaryElemInfo(fi).elem());

  const Real face_rho = _rho(boundary_face, time_arg);
  RealVectorValue density_times_velocity;
  for (const auto component : index_range(_vel))
    density_times_velocity(component) =
        boundary_normal_multiplier * face_rho * boundaryVelocityValue(fi, component, time_arg);

  return density_times_velocity * fi->normal();
}

Real
RhieChowMassFlux::boundaryVolumetricFluxTarget(const FaceInfo * fi,
                                               const Moose::StateArg & time_arg) const
{
  mooseAssert(fi && !_vel[0]->isInternalFace(*fi),
              "boundaryVolumetricFluxTarget should only be called on boundary faces.");

  const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);

  RealVectorValue face_velocity;
  for (const auto component : index_range(_vel))
    face_velocity(component) =
        boundary_normal_multiplier * boundaryVelocityValue(fi, component, time_arg);

  return face_velocity * fi->normal();
}

Real
RhieChowMassFlux::computeFaceNormalAinv(const RealVectorValue & face_ainv,
                                        const RealVectorValue & face_normal) const
{
  Real normal_ainv = 0.0;
  for (const auto dim_i : make_range(_dim))
    normal_ainv += face_ainv(dim_i) * face_normal(dim_i) * face_normal(dim_i);

  return normal_ainv;
}

Real
RhieChowMassFlux::boundaryNormalAinv(const FaceInfo * fi) const
{
  mooseAssert(fi && !_vel[0]->isInternalFace(*fi),
              "boundaryNormalAinv should only be called on boundary faces.");

  const auto & face_ainv = libmesh_map_find(_Ainv, fi->id());
  return computeFaceNormalAinv(face_ainv, fi->normal());
}

bool
RhieChowMassFlux::isAdjustablePressureBoundaryFace(const FaceInfo * fi) const
{
  if (!fi || _vel[0]->isInternalFace(*fi) || fi->boundaryIDs().empty())
    return false;

  if (const auto * bc_pointer = pressureBoundaryCondition(fi))
    return dynamic_cast<const LinearFVPressureFluxBC *>(bc_pointer) ||
           dynamic_cast<const LinearFVPressureSymmetryBC *>(bc_pointer);

  return false;
}

Real
RhieChowMassFlux::pressurePredictorFluxAdjustment(const FaceInfo * fi) const
{
  if (!fi)
    return 0.0;

  const auto it = _pressure_predictor_flux_adjustment.find(fi->id());
  return it == _pressure_predictor_flux_adjustment.end() ? 0.0 : it->second;
}

Real
RhieChowMassFlux::pressureBoundaryTargetFlux(const FaceInfo * fi,
                                             const Moose::StateArg & time_arg) const
{
  return boundaryMassFluxTarget(fi, time_arg);
}

Real
RhieChowMassFlux::pressureBoundaryNormalAinv(const FaceInfo * fi) const
{
  return boundaryNormalAinv(fi);
}

LinearFVBoundaryCondition *
RhieChowMassFlux::pressureBoundaryCondition(const FaceInfo * fi) const
{
  return fi ? _p->getBoundaryCondition(*fi) : nullptr;
}

Real
RhieChowMassFlux::pressureDiffusionFaceArea(const FaceInfo * /*fi*/) const
{
  return 1.0;
}

bool
RhieChowMassFlux::useConstrainedBoundaryPredictorState(const FaceInfo * fi) const
{
  if (!fi || _vel[0]->isInternalFace(*fi))
    return false;

  // A LinearFVPressureFluxBC already enforces the prescribed boundary mass flux, so reconstructing
  // the predictor state from the Dirichlet velocity there would impose it twice. Those faces take
  // the one-term boundary expansion instead.
  if (dynamic_cast<const LinearFVPressureFluxBC *>(pressureBoundaryCondition(fi)))
    return false;

  return _vel[0]->isDirichletBoundaryFace(*fi);
}

bool
RhieChowMassFlux::isFlowFace(const FaceInfo & fi) const
{
  return hasBlocks(fi.elemPtr()->subdomain_id()) ||
         (fi.neighborPtr() && hasBlocks(fi.neighborPtr()->subdomain_id()));
}

const ElemInfo &
RhieChowMassFlux::boundaryElemInfo(const FaceInfo * fi) const
{
  mooseAssert(fi && !_vel[0]->isInternalFace(*fi),
              "boundaryElemInfo should only be called on boundary faces.");
  return hasBlocks(fi->elemPtr()->subdomain_id()) ? *fi->elemInfo() : *fi->neighborInfo();
}

Real
RhieChowMassFlux::boundaryNormalMultiplier(const FaceInfo * fi) const
{
  mooseAssert(fi && !_vel[0]->isInternalFace(*fi),
              "boundaryNormalMultiplier should only be called on boundary faces.");
  return hasBlocks(fi->elemPtr()->subdomain_id()) ? 1.0 : -1.0;
}

Moose::FaceArg
RhieChowMassFlux::makeCenteredFaceArg(const FaceInfo * fi,
                                      const Elem * const face_side,
                                      const Moose::StateArg * limiter_state) const
{
  return Moose::FaceArg{
      fi, Moose::FV::LimiterType::CentralDifference, true, false, face_side, limiter_state};
}

Real
RhieChowMassFlux::faceMeasure(const FaceInfo & fi) const
{
  return fi.faceArea() * fi.faceCoord();
}

void
RhieChowMassFlux::updatePressureBoundaryNormalGradients(const bool apply_pressure_flux_adjustment)
{
  if (!_pressure_predictor_face_state_valid)
    updatePressurePredictorFaceState();

  const auto time_arg = Moose::currentState();
  Real max_boundary_normal_ainv = 0.0;

  for (const auto * fi : _flow_face_info)
  {
    _pressure_predictor_flux_adjustment[fi->id()] = 0.0;
    // Rebuild the live predictor source in the native pressure-correction
    // space so constrained pressure BCs see the same total phiHbyA contract
    // that the interior pressure equation sees. The diffuse-interface path stores
    // _pressure_predictor_base_flux with the internal sign, while _phig_flux is
    // the explicit body-force contribution, so internal phiHbyA is base - phig.
    _pressure_predictor_flux[fi->id()] =
        _pressure_predictor_base_flux[fi->id()] - _phig_flux[fi->id()];
    if (!_vel[0]->isInternalFace(*fi))
    {
      _pressure_boundary_normal_gradient[fi->id()] = 0.0;
      max_boundary_normal_ainv =
          std::max(max_boundary_normal_ainv, std::abs(pressureBoundaryNormalAinv(fi)));
    }
  }

  // Treat a boundary pressure coefficient as degenerate only relative to the scale of the
  // current boundary operator. An absolute libMesh::TOLERANCE cutoff is too aggressive for
  // high-density-ratio outlet faces, where valid normal Ainv values can be smaller than 1e-6.
  const Real degenerate_normal_ainv_tol =
      std::max(std::numeric_limits<Real>::min(),
               std::numeric_limits<Real>::epsilon() * std::max(1.0, max_boundary_normal_ainv));

  // Shift the predictor source on adjustable fixed-flux patches so the integrated predictor
  // boundary flux matches the boundary mass-flux target before the pressure solve.
  if (apply_pressure_flux_adjustment)
  {
    Real integrated_source_imbalance = 0.0;
    Real adjustable_measure = 0.0;

    for (const auto * fi : _flow_face_info)
      if (isAdjustablePressureBoundaryFace(fi))
      {
        const Real face_measure = faceMeasure(*fi);
        const Real required_pressure_flux =
            _pressure_predictor_flux[fi->id()] + pressureBoundaryTargetFlux(fi, time_arg);
        integrated_source_imbalance += pressureCorrectionFluxIsIntegrated()
                                           ? required_pressure_flux
                                           : required_pressure_flux * face_measure;
        adjustable_measure += face_measure;
      }

    if (adjustable_measure > libMesh::TOLERANCE)
    {
      const Real uniform_source_adjustment_density =
          -integrated_source_imbalance / adjustable_measure;
      for (const auto * fi : _flow_face_info)
        if (isAdjustablePressureBoundaryFace(fi))
        {
          const Real face_measure = faceMeasure(*fi);
          const Real source_adjustment = pressureCorrectionFluxIsIntegrated()
                                             ? uniform_source_adjustment_density * face_measure
                                             : uniform_source_adjustment_density;
          _pressure_predictor_flux_adjustment[fi->id()] = source_adjustment;
          _pressure_predictor_flux[fi->id()] += source_adjustment;
        }
    }
  }

  for (const auto * fi : _flow_face_info)
  {
    if (_vel[0]->isInternalFace(*fi))
      continue;

    const Real normal_ainv = pressureBoundaryNormalAinv(fi);
    if (!std::isfinite(normal_ainv) || std::abs(normal_ainv) <= degenerate_normal_ainv_tol)
    {
      _pressure_boundary_normal_gradient[fi->id()] = 0.0;
      continue;
    }

    const Real required_pressure_flux =
        _pressure_predictor_flux[fi->id()] + pressureBoundaryTargetFlux(fi, time_arg);
    const Real face_measure = faceMeasure(*fi);
    if (pressureCorrectionFluxIsIntegrated() && face_measure <= libMesh::TOLERANCE)
    {
      _pressure_boundary_normal_gradient[fi->id()] = 0.0;
      continue;
    }

    const Real required_pressure_flux_density = pressureCorrectionFluxIsIntegrated()
                                                    ? required_pressure_flux / face_measure
                                                    : required_pressure_flux;
    _pressure_boundary_normal_gradient[fi->id()] = -required_pressure_flux_density / normal_ainv;
  }

  _pressure_boundary_normal_gradient_valid = true;
  ++_face_mass_flux_generation;
}

void
RhieChowMassFlux::initCouplingField()
{
  // We loop through the faces and populate the coupling fields (face H/A and 1/H)
  // with 0s for now. Pressure corrector solves will always come after the
  // momentum source so we expect these fields to change before the actual solve.
  for (auto & fi : _fe_problem.mesh().faceInfo())
  {
    _Ainv[fi->id()];
    _HbyA_flux[fi->id()];
    _phiHbyA_flux[fi->id()];
    _phig_flux[fi->id()];
  }
}

RealVectorValue
RhieChowMassFlux::interpolateMomentumVectorReadersToFace(
    const FaceInfo * fi, const std::vector<PetscVectorReader> & readers) const
{
  using namespace Moose::FV;

  RealVectorValue face_value;

  if (!fi || readers.size() < _dim)
    return face_value;

  if (_vel[0]->isInternalFace(*fi))
  {
    const auto & elem_info = *fi->elemInfo();
    const auto & neighbor_info = *fi->neighborInfo();

    for (const auto dim_i : make_range(_dim))
    {
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      const auto neighbor_dof =
          neighbor_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      face_value(dim_i) =
          linearInterpolation(readers[dim_i](elem_dof), readers[dim_i](neighbor_dof), *fi, true);
    }

    return face_value;
  }

  const ElemInfo & elem_info = boundaryElemInfo(fi);
  for (const auto dim_i : make_range(_dim))
  {
    const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
    face_value(dim_i) = readers[dim_i](elem_dof);
  }

  return face_value;
}

RhieChowMassFlux::FaceMomentumPredictorState
RhieChowMassFlux::buildMomentumPredictorFaceState(
    const FaceInfo * fi,
    const std::vector<PetscVectorReader> & hbya_reader,
    const std::vector<PetscVectorReader> * ainv_reader,
    const Moose::StateArg & time_arg,
    const bool include_boundary_body_force_correction,
    const bool use_boundary_velocity_values) const
{
  using namespace Moose::FV;

  FaceMomentumPredictorState state;
  if (!fi)
    return state;

  mooseAssert(hbya_reader.size() >= _dim, "HbyA readers must cover every momentum component.");
  mooseAssert(!ainv_reader || ainv_reader->size() >= _dim,
              "Ainv readers must cover every momentum component.");

  if (_vel[0]->isInternalFace(*fi))
  {
    const auto & elem_info = *fi->elemInfo();
    const auto & neighbor_info = *fi->neighborInfo();
    const Real elem_rho = _rho(makeElemArg(fi->elemPtr()), time_arg);
    const Real neighbor_rho = _rho(makeElemArg(fi->neighborPtr()), time_arg);

    interpolate(
        Moose::FV::InterpMethod::Average, state.face_rho, elem_rho, neighbor_rho, *fi, true);
    if (ainv_reader)
      state.face_pressure_ainv = interpolateMomentumVectorReadersToFace(fi, *ainv_reader);

    for (const auto dim_i : make_range(_dim))
    {
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      const auto neighbor_dof =
          neighbor_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      interpolate(Moose::FV::InterpMethod::Average,
                  state.face_hbya(dim_i),
                  hbya_reader[dim_i](elem_dof),
                  hbya_reader[dim_i](neighbor_dof),
                  *fi,
                  true);
      interpolate(Moose::FV::InterpMethod::Average,
                  state.density_weighted_face_hbya(dim_i),
                  elem_rho * hbya_reader[dim_i](elem_dof),
                  neighbor_rho * hbya_reader[dim_i](neighbor_dof),
                  *fi,
                  true);

      if (ainv_reader)
      {
        interpolate(_pressure_diffusion_interp_method,
                    state.face_density_weighted_ainv(dim_i),
                    elem_rho * (*ainv_reader)[dim_i](elem_dof),
                    neighbor_rho * (*ainv_reader)[dim_i](neighbor_dof),
                    *fi,
                    true);
      }
    }

    return state;
  }

  const Real boundary_normal_multiplier = boundaryNormalMultiplier(fi);
  const ElemInfo & elem_info = boundaryElemInfo(fi);
  const Moose::FaceArg boundary_face = makeCenteredFaceArg(fi, elem_info.elem());
  const bool use_constrained_boundary_state = useConstrainedBoundaryPredictorState(fi);

  if (use_constrained_boundary_state)
  {
    state.face_rho = _rho(boundary_face, Moose::currentState());
    for (const auto dim_i : make_range(_dim))
    {
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      Real boundary_value = std::numeric_limits<Real>::quiet_NaN();
      if (use_boundary_velocity_values)
        boundary_value = boundaryVelocityValue(fi, dim_i, time_arg);
      else if (auto * bc_pointer = velocityBoundaryCondition(fi, dim_i))
        boundary_value = bc_pointer->computeBoundaryValue();

      state.face_hbya(dim_i) =
          std::isfinite(boundary_value)
              ? -boundary_value
              : -MetaPhysicL::raw_value((*_vel[dim_i])(boundary_face, Moose::currentState()));

      if (include_boundary_body_force_correction && !_split_momentum_predictor_operator &&
          !_body_force_kernel_names.empty())
        for (const auto & force_kernel : _body_force_kernels[dim_i])
        {
          force_kernel->setCurrentElemInfo(&elem_info);
          state.face_hbya(dim_i) -=
              force_kernel->computeRightHandSideContribution() * (*ainv_reader)[dim_i](elem_dof) /
              (elem_info.volume() * elem_info.coordFactor()); // zero-term expansion
        }

      state.face_hbya(dim_i) *= boundary_normal_multiplier;
      state.density_weighted_face_hbya(dim_i) = state.face_rho * state.face_hbya(dim_i);
    }
  }
  else
  {
    state.face_rho = _rho(makeElemArg(elem_info.elem()), time_arg);
    for (const auto dim_i : make_range(_dim))
    {
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      state.face_hbya(dim_i) = boundary_normal_multiplier * hbya_reader[dim_i](elem_dof);
      state.density_weighted_face_hbya(dim_i) = state.face_rho * state.face_hbya(dim_i);
    }
  }

  if (ainv_reader)
  {
    state.face_pressure_ainv = interpolateMomentumVectorReadersToFace(fi, *ainv_reader);
    const Real elem_rho = _rho(makeElemArg(elem_info.elem()), time_arg);
    for (const auto dim_i : make_range(_dim))
    {
      const auto elem_dof = elem_info.dofIndices()[_global_momentum_system_numbers[dim_i]][0];
      state.face_density_weighted_ainv(dim_i) = elem_rho * (*ainv_reader)[dim_i](elem_dof);
    }
  }

  return state;
}

void
RhieChowMassFlux::populateCouplingFunctors(
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

  // We loop through the faces and populate the coupling fields (face H/A and 1/A)
  for (auto & fi : _flow_face_info)
  {
    const auto state =
        buildMomentumPredictorFaceState(fi,
                                        hbya_reader,
                                        &ainv_reader,
                                        time_arg,
                                        /* include_boundary_body_force_correction = */
                                        true,
                                        /* use_boundary_velocity_values = */ false);
    _Ainv[fi->id()] = state.face_density_weighted_ainv;
    _pressure_Ainv[fi->id()] = state.face_pressure_ainv;
    // Lastly, we populate the face flux resulted by H/A
    _HbyA_flux[fi->id()] = state.density_weighted_face_hbya * fi->normal();
  }
}

void
RhieChowMassFlux::computeHbyA(bool verbose)
{
  if (verbose)
  {
    _console << "************************************" << std::endl;
    _console << "Computing HbyA" << std::endl;
    _console << "************************************" << std::endl;
  }
  mooseAssert(_momentum_implicit_systems.size() && _momentum_implicit_systems[0],
              "The momentum system shall be linked before calling this function!");

  mooseAssert(!_grad_p_current.empty(),
              "A coupling pressure-gradient snapshot must exist before computing H/A.");
  mooseAssert(_grad_p_current.size() == _momentum_implicit_systems.size(),
              "The coupling pressure-gradient snapshot must have one component per momentum "
              "system.");

#ifndef NDEBUG
  for (const auto system_i : index_range(_momentum_implicit_systems))
  {
    const auto & snapshot = *_grad_p_current[system_i];
    const auto & momentum_layout = *_momentum_implicit_systems[system_i]->current_local_solution;
    mooseAssert(snapshot.size() == momentum_layout.size() &&
                    snapshot.local_size() == momentum_layout.local_size() &&
                    snapshot.first_local_index() == momentum_layout.first_local_index() &&
                    snapshot.last_local_index() == momentum_layout.last_local_index(),
                "The coupling pressure-gradient snapshot layout must match its momentum system.");
  }
#endif

  const auto & coupling_pressure_gradient = _grad_p_current;
  _pressure_equation_flux_valid = false;
  _pressure_boundary_normal_gradient_valid = false;
  _pressure_predictor_face_state_valid = false;

  const bool use_cached_predictor_operator = canUseCachedMomentumPredictorOperator();
  const bool split_predictor_operator = _split_momentum_predictor_operator;

  if (verbose && use_cached_predictor_operator)
    _console << "Using cached momentum predictor operator for HbyA/Ainv." << std::endl;

  _HbyA_raw.clear();
  _Ainv_raw.clear();

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

    if (verbose)
    {
      _console << "Velocity solution in H(u)" << std::endl;
      solution.print();
    }

    if (use_cached_predictor_operator)
    {
      _Ainv_raw.push_back(current_local_solution.zero_clone());
      _HbyA_raw.push_back(current_local_solution.zero_clone());
      const auto * const predictor = cachedMomentumPredictorOperator(system_i);
      mooseAssert(predictor && predictor->complete(),
                  "A complete momentum predictor operator is required before computing HbyA.");

      computePredictorOperatorBase(
          system_i, *(_HbyA_raw.back()), *(_Ainv_raw.back()), predictor->rhs_raw.get());

      auto h_correction = current_local_solution.zero_clone();
      auto * h_correction_petsc = dynamic_cast<PetscVector<Number> *>(h_correction.get());
      mooseAssert(h_correction_petsc,
                  "The vectors used in the RhieChowMassFlux objects need to be convertable "
                  "to PetscVectors!");
      h_correction_petsc->pointwise_mult(*predictor->h_diagonal_correction_raw,
                                         current_local_solution);
      h_correction_petsc->close();
      _HbyA_raw.back()->add(1.0, *h_correction);

      *(_Ainv_raw.back()) = *predictor->pressure_coupling_diagonal_raw;
      _Ainv_raw.back()->close();
    }
    else
    {
      _Ainv_raw.push_back(current_local_solution.zero_clone());
      _HbyA_raw.push_back(current_local_solution.zero_clone());
      computePredictorOperatorBase(system_i, *(_HbyA_raw.back()), *(_Ainv_raw.back()), nullptr);
    }

    NumericVector<Number> & Ainv = *(_Ainv_raw.back());
    NumericVector<Number> & HbyA = *(_HbyA_raw.back());

    // We create a working vector to ease some of the operations.
    auto working_vector = momentum_system->current_local_solution->zero_clone();
    PetscVector<Number> * working_vector_petsc =
        dynamic_cast<PetscVector<Number> *>(working_vector.get());
    mooseAssert(working_vector_petsc,
                "The vectors used in the RhieChowMassFlux objects need to be convertable "
                "to PetscVectors!");

    if (!split_predictor_operator)
    {
      // The non-split predictor operator includes explicit pressure/body-force
      // terms in the assembled RHS, so HbyA has to remove pressure and then
      // restore explicit body-force sources.
      // The coupling gradient is the snapshot taken before this predictor was assembled,
      // which is the gradient the momentum pressure source actually used.
      working_vector_petsc->pointwise_mult(*coupling_pressure_gradient[system_i], *_cell_volumes);
      HbyA.add(-1.0, *working_vector_petsc);

      if (!_body_force_kernels.empty())
      {
        auto explicit_body_force_rhs = current_local_solution.zero_clone();
        for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
          if (hasBlocks(elem_info->subdomain_id()))
          {
            const auto elem_dof =
                elem_info->dofIndices()[_global_momentum_system_numbers[system_i]][0];

            Real rhs_contribution = 0.0;
            for (const auto & force_kernel : _body_force_kernels[system_i])
            {
              force_kernel->setCurrentElemInfo(elem_info);
              rhs_contribution += force_kernel->computeRightHandSideContribution();
            }

            explicit_body_force_rhs->set(elem_dof, rhs_contribution);
          }

        explicit_body_force_rhs->close();
        HbyA.add(1.0, *explicit_body_force_rhs);
      }
    }

    if (verbose)
    {
      _console << "total RHS" << std::endl;
      rhs.print();
      _console << "pressure RHS" << std::endl;
      coupling_pressure_gradient[system_i]->print();
      _console << " H(u)-rhs-relaxation_source" << std::endl;
      if (!split_predictor_operator)
      {
        _console << "pressure RHS" << std::endl;
        coupling_pressure_gradient[system_i]->print();
        _console << " predictor_base - pressure + explicit_body_force" << std::endl;
      }
      else
        _console << " predictor_base (split operator)" << std::endl;
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
      working_vector_petsc->pointwise_mult(*Ainv_full, *coupling_pressure_gradient[system_i]);
      working_vector_petsc->pointwise_mult(*working_vector_petsc, *_cell_volumes);
      HbyA.add(-1.0, *working_vector_petsc);

      // Correct Ainv
      Ainv = *Ainv_full_old;
    }

    Ainv.pointwise_mult(Ainv, *_cell_volumes);

    if (verbose)
    {
      _console << " 1/A" << std::endl;
      Ainv.print();
    }
  }

  // We fill the 1/A and H/A functors
  populateCouplingFunctors(_HbyA_raw, _Ainv_raw);

  if (verbose)
  {
    _console << "************************************" << std::endl;
    _console << "DONE Computing HbyA " << std::endl;
    _console << "************************************" << std::endl;
  }
}

const std::vector<NumericVector<Number> *> &
RhieChowMassFlux::pressureGradientComponents() const
{
  return pressureGradientField().components();
}

bool
RhieChowMassFlux::usingReconstructedPressureGradientMethod() const
{
  return dynamic_cast<const FVReconstructedPressureGradient *>(&pressureGradientField().method());
}

const FVReconstructedPressureGradient &
RhieChowMassFlux::reconstructedGradientMethod() const
{
  mooseAssert(usingReconstructedPressureGradientMethod(),
              "reconstructedGradientMethod() should only be called when reconstructed "
              "pressure gradients are active.");

  return dynamic_cast<const FVReconstructedPressureGradient &>(pressureGradientField().method());
}

FVReconstructedPressureGradient &
RhieChowMassFlux::reconstructedGradientMethod()
{
  auto & method = _fe_problem.getFVGradientMethod(pressureGradientField().method().name(), _tid);
  mooseAssert(&method == &pressureGradientField().method(),
              "The registered writable gradient method must match the pressure-gradient reader.");
  return dynamic_cast<FVReconstructedPressureGradient &>(method);
}

void
RhieChowMassFlux::checkReconstructedPressureGradientCompatibility() const
{
  if (!_pressure_system || !_pressure_gradient_field || !_base_pressure_gradient_field ||
      !_p_diffusion_kernel)
    mooseError("RhieChowMassFlux '",
               name(),
               "' is missing a required pressure system, pressure-gradient field, base-gradient "
               "field, or pressure-correction diffusion kernel for reconstructed pressure "
               "gradients.");

  reconstructedGradientMethod().validateSetup(*this);

  if (_momentum_systems.size() != _dim || _momentum_implicit_systems.size() != _dim ||
      _global_momentum_system_numbers.size() != _dim)
    mooseError("RhieChowMassFlux '",
               name(),
               "' requires one linked momentum system per spatial component for reconstructed "
               "pressure gradients.");

  if (_pressure_system->nVariables() != 1)
    mooseError(
        "FVReconstructedPressureGradient assumes the pressure and momentum systems each have "
        "exactly one variable so their DOF indices can be used in vector operations. Pressure "
        "system '",
        _pressure_system->name(),
        "' has variables: ",
        Moose::stringify(_pressure_system->getVariableNames()));

  for (const auto system_i : index_range(_momentum_systems))
  {
    const auto * const momentum_system = _momentum_systems[system_i];
    if (!momentum_system || !_momentum_implicit_systems[system_i])
      mooseError("RhieChowMassFlux '",
                 name(),
                 "' has an invalid momentum-system layout for component ",
                 system_i,
                 ".");

    if (momentum_system->nVariables() != 1)
      mooseError("FVReconstructedPressureGradient assumes the pressure and momentum systems each "
                 "have exactly one variable so their DOF indices can be used in vector operations. "
                 "Momentum system '",
                 momentum_system->name(),
                 "' has variables: ",
                 Moose::stringify(momentum_system->getVariableNames()));
  }
}
