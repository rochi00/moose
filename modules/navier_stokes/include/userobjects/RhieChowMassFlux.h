//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RhieChowFaceFluxProvider.h"
#include "CellCenteredMapFunctor.h"
#include "FaceCenteredMapFunctor.h"
#include "VectorComponentFunctor.h"
#include "LinearFVElementalKernel.h"
#include "LinearFVAnisotropicDiffusion.h"
#include "PetscVectorReader.h"
#include <string>
#include <unordered_map>
#include <set>

#include "libmesh/petsc_vector.h"

class MooseMesh;
class ElemInfo;
class INSFVVelocityVariable;
class INSFVPressureVariable;
class LinearFVGradientReader;
class FVReconstructedPressureGradient;
class LinearFVPressureCorrectionDiffusion;
class LinearFVBoundaryCondition;
class LinearFVElementalKernel;
namespace libMesh
{
class Elem;
class MeshBase;
template <typename T>
class SparseMatrix;
}

/**
 * User object responsible for determining the face fluxes using the Rhie-Chow interpolation in a
 * segregated solver that uses the linear FV formulation.
 */
class RhieChowMassFlux : public RhieChowFaceFluxProvider, public NonADFunctorInterface
{
public:
  static InputParameters validParams();
  RhieChowMassFlux(const InputParameters & params);
  ~RhieChowMassFlux() override;

  /// Get the face velocity times density (used in advection terms)
  virtual Real getMassFlux(const FaceInfo & fi) const;

  /// Get the volumetric face flux (used in advection terms)
  virtual Real getVolumetricFaceFlux(const FaceInfo & fi) const;

  /// Get the registered pressure gradient field used by compatible momentum pressure kernels.
  const LinearFVGradientReader & pressureGradientField() const;

  /// Get the base pressure gradient field used to seed reconstructed gradient updates.
  const LinearFVGradientReader & basePressureGradientField() const;

  /// Get the base-pressure gradient component vectors used to seed reconstructed updates.
  const std::vector<NumericVector<Number> *> & basePressureGradientComponents() const;

  /// Get the momentum-layout H/A component vectors.
  const std::vector<std::unique_ptr<NumericVector<Number>>> & HbyAComponents() const;

  /// Get the momentum-layout 1/A component vectors.
  const std::vector<std::unique_ptr<NumericVector<Number>>> & AinvComponents() const;

  /// Variable number of the pressure variable reconstructed by this object.
  unsigned int pressureVariableNumber() const;

  /// Number of spatial dimensions of the coupled momentum/pressure systems.
  unsigned int dimension() const { return _dim; }

  /// Velocity variable for one spatial component.
  MooseLinearVariableFVReal & velocityVariable(unsigned int component) { return *_vel[component]; }

  /// Velocity variable for one spatial component.
  const MooseLinearVariableFVReal & velocityVariable(unsigned int component) const
  {
    return *_vel[component];
  }

  /// Momentum linear system for one spatial component.
  LinearSystem & momentumSystem(unsigned int component) { return *_momentum_systems[component]; }

  /// Momentum linear system for one spatial component.
  const LinearSystem & momentumSystem(unsigned int component) const
  {
    return *_momentum_systems[component];
  }

  /// Pressure linear system.
  const LinearSystem & pressureSystem() const { return *_pressure_system; }

  /**
   * Producer counter for the conservative face mass flux: incremented every time
   * computeFaceMassFlux() runs. Used to verify reconstruction consumes the face flux produced by
   * the current pressure corrector, not one left over from a preceding corrector.
   */
  dof_id_type faceMassFluxGeneration() const { return _face_mass_flux_generation; }

  virtual Real getVolumetricFaceFlux(const Moose::FV::InterpMethod m,
                                     const FaceInfo & fi,
                                     const Moose::StateArg & time,
                                     const THREAD_ID tid,
                                     bool subtract_mesh_velocity) const override;

  /// Initialize the container for face velocities
  void initFaceMassFlux();
  /// Initialize the coupling fields (HbyA and Ainv)
  void initCouplingField();
  /// Cache the pressure-equation face flux from the current solved pressure field.
  virtual void cachePressureEquationFlux();
  /// Maximum face-flux Courant number over active flow cells for the supplied timestep.
  Real maxCourant(const Real dt) const;
  /// Update the values of the face velocities in the containers
  virtual void computeFaceMassFlux();
  /**
   * Writes the pressure-corrected cell velocity. This is the seam the corrector goes through, so
   * that a provider which reconstructs that velocity differently is given the chance to; the
   * arithmetic for the ordinary case lives in updateCellVelocity.
   */
  virtual void computeCellVelocity();
  /// Accessor for the current cell Ainv state.
  Real cellAinvRaw(const unsigned int system_i, const dof_id_type dof) const;
  /// Begin streaming native linear FV assembly contributions into the cached split predictor.
  void beginFVSplitMomentumPredictorOperatorAssembly(const unsigned int system_i);
  /// Finalize the split predictor assembled from the momentum operator.
  void completeFVSplitMomentumPredictorOperatorAssembly(const unsigned int system_i,
                                                        const Real relaxation_parameter,
                                                        const bool enforce_diagonal_dominance);
  /// Populate cached split-predictor internals from the assembled momentum operator.
  void cacheFVSplitMomentumPredictorOperatorAssembly(const unsigned int system_i,
                                                     libMesh::SparseMatrix<Number> & matrix,
                                                     const NumericVector<Number> & rhs);
  /// Invalidate the cached assembled/relaxed momentum predictor operator.
  void clearMomentumPredictorOperatorCache();
  /// Add explicit forcing to the momentum predictor RHS after the base operator assembly.
  virtual void addMomentumPredictorExplicitForcing(const unsigned int system_i,
                                                   NumericVector<Number> & rhs) const;
  /// Whether the predictor operator excludes explicit pressure/body-force terms.
  bool splitMomentumPredictorOperator() const { return _split_momentum_predictor_operator; }
  /// Update boundary pressure gradients from the current predictor and boundary velocity state.
  void updatePressureBoundaryNormalGradients(const bool apply_pressure_flux_adjustment);
  /// Refresh boundary-face velocity values from the active FV velocity BC objects.
  void updateVelocityBoundaryState();

  /// Whether the registered pressure gradient field is produced by the reconstructed method.
  bool usingReconstructedPressureGradientMethod() const;

  /// Snapshot the fields read by the next momentum predictor.
  void prepareMomentumPredictor();

  /// Capture the lagged velocity gradients computed while assembling the momentum predictor.
  void finalizeMomentumPredictor();

  /// Capture the lagged velocity gradient for another corrector using the current predictor.
  void preparePISOCorrector();

  /**
   * Before pressure relaxation, form the conservative reconstructed candidate from the corrected
   * face flux and update the cell velocity from that unrelaxed candidate. No-op for ordinary
   * pressure-gradient methods.
   */
  void preparePressureRelaxation();

  /**
   * After pressure relaxation, refresh every registered pressure gradient and publish relaxed
   * reconstructed feedback for the next momentum predictor. Ordinary methods also update the cell
   * velocity here, preserving their existing relaxed-pressure behavior.
   */
  void finalizePressureCorrector();

  /// Compute every registered pressure gradient once, before the first momentum predictor of a
  /// SIMPLE iteration sequence has anything to read.
  void initPressureGradient();

  virtual void meshChanged() override;
  virtual void initialize() override;
  virtual void execute() override {}
  virtual void finalize() override {}
  virtual void initialSetup() override;

  /**
   * Prepare reconstructed-gradient state once per attempted time step. Accepted coupling feedback
   * is retained for the next time step and restored when an attempt is retried, while candidate
   * fields and sequencing counters are reset for every attempt.
   */
  virtual void timestepSetup() override;

  /**
   * Update the momentum system-related information
   * @param momentum_systems Pointers to the momentum systems which are solved for the momentum
   * vector components
   * @param pressure_system Reference to the pressure system
   * @param momentum_system_numbers The numbers of these systems
   */
  void linkMomentumPressureSystems(const std::vector<LinearSystem *> & momentum_systems,
                                   LinearSystem & pressure_system,
                                   const std::vector<unsigned int> & momentum_system_numbers);

  /**
   * Computes the inverse of the diagonal (1/A) of the system matrix plus the H/A components for the
   * pressure equation plus Rhie-Chow interpolation.
   */
  void computeHbyA(const bool verbose);

protected:
  /// Update cell velocity from the supplied momentum-coupling pressure gradient.
  template <typename GradientComponent>
  void updateCellVelocity(const std::vector<GradientComponent> & pressure_gradient);

  /// Get the registered pressure gradient component vectors.
  const std::vector<NumericVector<Number> *> & pressureGradientComponents() const;

  /// Access the reconstructed gradient method when it is selected.
  FVReconstructedPressureGradient & reconstructedGradientMethod();

  /// Access the reconstructed gradient method when it is selected.
  const FVReconstructedPressureGradient & reconstructedGradientMethod() const;

  /// Check the single-variable system layout assumed by reconstructed pressure-gradient vector ops.
  void checkReconstructedPressureGradientCompatibility() const;
  struct MomentumPredictorOperator;
  struct FaceMomentumPredictorState
  {
    Real face_rho = 0.0;
    RealVectorValue face_hbya;
    RealVectorValue density_weighted_face_hbya;
    RealVectorValue face_density_weighted_ainv;
    RealVectorValue face_pressure_ainv;
  };

  /// Predictor-side face flux entering the pressure equation. Defaults to the HbyA contribution.
  Real pressurePredictorFlux(const FaceInfo * fi) const;
  /// Populate the explicit face state used by the pressure corrector.
  void updatePressurePredictorFaceState();
  /// Optional per-face adjustment added to the predictor flux for reference-pressure compatibility.
  Real pressurePredictorFluxAdjustment(const FaceInfo * fi) const;
  /// Target boundary flux for the active pressure-correction space.
  virtual Real pressureBoundaryTargetFlux(const FaceInfo * fi,
                                          const Moose::StateArg & time_arg) const;
  /// Normal diffusion coefficient for the active pressure-correction space.
  virtual Real pressureBoundaryNormalAinv(const FaceInfo * fi) const;
  /// Fetch and initialize the pressure boundary condition attached to a boundary face.
  LinearFVBoundaryCondition * pressureBoundaryCondition(const FaceInfo * fi) const;
  /// Whether pressure-correction face fluxes include face area.
  virtual bool pressureCorrectionFluxIsIntegrated() const { return false; }
  /// Whether a boundary face should use constrained boundary values for the predictor state.
  virtual bool useConstrainedBoundaryPredictorState(const FaceInfo * fi) const;
  /// Whether a face touches the active flow block set.
  bool isFlowFace(const FaceInfo & fi) const;
  /// Flow-side element information for a boundary face.
  const ElemInfo & boundaryElemInfo(const FaceInfo * fi) const;
  /// Orientation multiplier from flow-side boundary values to face normal convention.
  Real boundaryNormalMultiplier(const FaceInfo * fi) const;
  /// Build a central-difference face argument on the requested side.
  Moose::FaceArg makeCenteredFaceArg(const FaceInfo * fi,
                                     const Elem * face_side = nullptr,
                                     const Moose::StateArg * limiter_state = nullptr) const;
  /// Geometric face measure used to convert between flux density and integrated flux.
  Real faceMeasure(const FaceInfo & fi) const;
  /// Face area supplied to the pressure diffusion kernel when reconstructing pressure fluxes.
  virtual Real pressureDiffusionFaceArea(const FaceInfo * fi) const;
  /// Evaluate or fetch the boundary face value for a velocity component.
  Real boundaryVelocityValue(const FaceInfo * fi,
                             const unsigned int component,
                             const Moose::StateArg & time_arg) const;
  /// Fetch and initialize the velocity boundary condition for one component.
  LinearFVBoundaryCondition * velocityBoundaryCondition(const FaceInfo * fi,
                                                        const unsigned int component) const;
  /// Compute the target physical boundary mass flux rho U_b.n at a face.
  Real boundaryMassFluxTarget(const FaceInfo * fi, const Moose::StateArg & time_arg) const;
  /// Compute the target volumetric boundary flux U_b.n at a face.
  Real boundaryVolumetricFluxTarget(const FaceInfo * fi, const Moose::StateArg & time_arg) const;
  /// Compute the normal component of a face diagonal Ainv vector.
  Real computeFaceNormalAinv(const RealVectorValue & face_ainv,
                             const RealVectorValue & face_normal) const;
  /// Compute the normal component of the face diffusion coefficient used by pressure BCs.
  Real boundaryNormalAinv(const FaceInfo * fi) const;
  /// Determine whether a boundary face belongs to an adjustable fixed-flux pressure patch.
  bool isAdjustablePressureBoundaryFace(const FaceInfo * fi) const;
  /// Compute the face flux contributed by the solved pressure equation in the
  /// native pressure-correction space used by this Rhie-Chow object.
  Real computeDiscretePressureFaceFlux(const FaceInfo * fi) const;
  /// Optional extension point after boundary velocity state has been refreshed.
  virtual void postUpdateVelocityBoundaryState() {}

  /// Accessor for the cached set of flow faces used by this Rhie-Chow object.
  const std::vector<const FaceInfo *> & flowFaceInfo() const { return _flow_face_info; }

  /// Compute the cell volumes on the mesh
  void setupMeshInformation();

  /// Populate the face values of the H/A and 1/A fields
  void
  populateCouplingFunctors(const std::vector<std::unique_ptr<NumericVector<Number>>> & raw_hbya,
                           const std::vector<std::unique_ptr<NumericVector<Number>>> & raw_Ainv);

  /// Build the face predictor state shared by stock and diffuse-interface pressure coupling.
  FaceMomentumPredictorState
  buildMomentumPredictorFaceState(const FaceInfo * fi,
                                  const std::vector<PetscVectorReader> & hbya_reader,
                                  const std::vector<PetscVectorReader> * ainv_reader,
                                  const Moose::StateArg & time_arg,
                                  bool include_boundary_body_force_correction,
                                  bool use_boundary_velocity_values) const;
  /// Interpolate per-component momentum-system cell readers to a face vector.
  RealVectorValue
  interpolateMomentumVectorReadersToFace(const FaceInfo * fi,
                                         const std::vector<PetscVectorReader> & readers) const;

  /// Build the base predictor operator vectors M*u - A*u - rhs and the relaxed diagonal.
  void computePredictorOperatorBase(const unsigned int system_i,
                                    NumericVector<Number> & base_raw,
                                    NumericVector<Number> & diagonal_raw,
                                    const NumericVector<Number> * rhs_override = nullptr) const;

  /// Refresh the aggregate validity flag for the cached predictor operator.
  void updateCachedMomentumPredictorOperatorValidity();

  /// Whether the cached predictor-operator state is complete and can be consumed.
  bool canUseCachedMomentumPredictorOperator();

  /// Finalize all component operators as one vector momentum-predictor object.
  void finalizeCachedMomentumPredictorOperators();

  /// Access the cached momentum-predictor operator for a component.
  const MomentumPredictorOperator *
  cachedMomentumPredictorOperator(const unsigned int system_i) const;

  /// Access the cached fvMatrix::D() diagonal used by pressure projection for a component.
  const NumericVector<Number> *
  cachedMomentumPredictorDiagonalRaw(const unsigned int system_i) const;

  /**
   * Check the block consistency between the passed in \p var and us
   */
  template <typename VarType>
  void checkBlocks(const VarType & var) const;

  virtual bool supportMeshVelocity() const override { return false; }

  /// The \p MooseMesh that this user object operates on
  const MooseMesh & _moose_mesh;

  /// The \p libMesh mesh that this object acts on
  const libMesh::MeshBase & _mesh;

  /// The dimension of the mesh, e.g. 3 for hexes and tets, 2 for quads and tris
  const unsigned int _dim;

  /// The thread 0 copy of the pressure variable
  const MooseLinearVariableFVReal * const _p;

  /// Thread 0 copies of the velocity variables; non-const to allow writing
  /// cell velocities in updateCellVelocity().
  std::vector<MooseLinearVariableFVReal *> _vel;

  /// Pointer to the pressure diffusion term in the pressure Poisson equation
  LinearFVAnisotropicDiffusion * _p_diffusion_kernel;

  /**
   * A map functor from faces to $HbyA_{ij} = (A_{offdiag}*\mathrm{(predicted~velocity)} -
   * \mathrm{Source})_{ij}/A_{ij}$. So this contains the off-diagonal part of the system matrix
   * multiplied by the predicted velocity minus the source terms from the right hand side of the
   * linearized momentum predictor step.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _HbyA_flux;

  /**
   * Explicit predictor-side face flux entering the pressure equation in the native
   * pressure-correction sign convention. The accepted transport-face flux is therefore formed as
   * -_phiHbyA_flux + _pressure_equation_flux.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _phiHbyA_flux;

  /**
   * Predictor-side face flux in the same native pressure-correction space consumed by the
   * pressure equation and pressure boundary-condition updates. For stock RC this is a normal
   * flux density; diffuse/conservative paths publish an area-integrated phi. It is
   * still a pressure-correction-space predictor quantity, not the accepted transport phi.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _pressure_predictor_flux;

  /**
   * Explicit face-force contribution entering the pressure corrector.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _phig_flux;

  /**
   * We hold on to the cell-based HbyA vectors so that we can easily reconstruct the
   * cell velocities as well.
   */
  std::vector<std::unique_ptr<NumericVector<Number>>> _HbyA_raw;

  /**
   * A map functor from faces to $(1/A)_f$. Where $A_i$ is the diagonal of the system matrix
   * for the momentum equation.
   */
  FaceCenteredMapFunctor<RealVectorValue, std::unordered_map<dof_id_type, RealVectorValue>> _Ainv;

  /**
   * Face diffusion tensor in the native pressure-correction space. For stock RC this matches the
   * density-weighted face Ainv, while diffuse/conservative paths may publish raw face rAUf here.
   */
  FaceCenteredMapFunctor<RealVectorValue, std::unordered_map<dof_id_type, RealVectorValue>>
      _pressure_Ainv;

  /**
   * We hold on to the cell-based 1/A vectors so that we can easily reconstruct the
   * cell velocities as well.
   */
  std::vector<std::unique_ptr<NumericVector<Number>>> _Ainv_raw;

  std::unique_ptr<NumericVector<Number>> _A_avg;

  /**
   * A map functor from faces to mass fluxes which are used in the advection terms.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> & _face_mass_flux;

  /// Coupling pressure gradient captured before the current momentum predictor is assembled.
  std::vector<std::unique_ptr<NumericVector<Number>>> _grad_p_current;
  /**
   * Exact face flux contribution from the solved pressure equation, in the same sign convention
   * used in face_mass_flux = -predictor_flux + pressure_equation_flux.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _pressure_equation_flux;

  /**
   * Cached boundary-normal pressure gradient used by fixed-flux pressure boundary updates against
   * the current predictor flux.
   */
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>>
      _pressure_boundary_normal_gradient;

  /// Indicates whether _pressure_equation_flux matches the latest solved pressure field.
  bool _pressure_equation_flux_valid = false;

  /// Per-face predictor source-flux adjustments used for the local adjustPhi analogue.
  std::unordered_map<dof_id_type, Real> _pressure_predictor_flux_adjustment;

  /// Unadjusted predictor-side face flux prior to the local adjustPhi analogue, stored in the
  /// internal pressure-correction convention. The accepted transport flux uses the negated value.
  std::unordered_map<dof_id_type, Real> _pressure_predictor_base_flux;

  /// Indicates whether _pressure_boundary_normal_gradient matches the latest predictor state.
  bool _pressure_boundary_normal_gradient_valid = false;

  /// Indicates whether _phiHbyA_flux/_phig_flux match the latest predictor state.
  bool _pressure_predictor_face_state_valid = false;

  /// Boundary-face velocity values refreshed from the current FV velocity BCs.
  std::vector<std::unordered_map<dof_id_type, Real>> _boundary_velocity_face_values;

  /// Indicates whether _boundary_velocity_face_values matches the latest cell-centered velocity.
  bool _velocity_boundary_state_valid = false;

  /// Whether to consume a cached assembled predictor operator instead of rebuilding HbyA live.
  const bool _use_cached_momentum_predictor_operator;

  /// Whether the momentum predictor operator already excludes explicit pressure/body-force terms.
  const bool _split_momentum_predictor_operator;

  /// Predictor operators for each momentum component.
  std::vector<std::unique_ptr<MomentumPredictorOperator>> _momentum_predictor_operators;

  /// Indicates whether the cached predictor operator has been populated for every component.
  bool _cached_predictor_operator_valid = false;

  /// Indicates whether all streamed component operators have been vector-finalized.
  bool _cached_predictor_operator_finalized = false;

  /// Pointer to the body force terms
  std::vector<std::vector<LinearFVElementalKernel *>> _body_force_kernels;
  /// Vector of body force term names
  std::vector<std::vector<std::string>> _body_force_kernel_names;

  /**
   * Producer counter for the conservative face mass flux: incremented every time
   * computeFaceMassFlux() runs. Exposed read-only through faceMassFluxGeneration() so the
   * reconstructed gradient can verify it is reconstructing from the face flux produced by the
   * current pressure corrector, not one left over from a preceding corrector.
   */
  dof_id_type _face_mass_flux_generation = 0;

  /**
   * Functor describing the density of the fluid
   */
  const Moose::Functor<Real> & _rho;

  /// Pointers to the linear system(s) in moose corresponding to the momentum equation(s)
  std::vector<LinearSystem *> _momentum_systems;

  /// Numbers of the momentum system(s)
  std::vector<unsigned int> _momentum_system_numbers;

  /// Global numbers of the momentum system(s)
  std::vector<unsigned int> _global_momentum_system_numbers;

  /// Pointers to the momentum equation implicit system(s) from libmesh
  std::vector<libMesh::LinearImplicitSystem *> _momentum_implicit_systems;

  /// Pointer to the pressure system
  LinearSystem * _pressure_system;

  /// Registered pressure gradient field used by Rhie-Chow and compatible momentum pressure kernels.
  const LinearFVGradientReader * _pressure_gradient_field;

  /// Base pressure gradient field used to seed reconstructed gradient updates.
  const LinearFVGradientReader * _base_pressure_gradient_field;

  /// Global number of the pressure system
  unsigned int _global_pressure_system_number;

  /// We will hold a vector of cell volumes to make sure we can do volume corrections rapidly
  std::unique_ptr<NumericVector<Number>> _cell_volumes;

  /// Enumerator for the method used for pressure projection
  const MooseEnum _pressure_projection_method;

  /**
   * Whether this provider maintains its own volumetric face flux rather than recovering one by
   * dividing the mass flux by a single face density. Only such a provider may interpolate the
   * density-velocity product, because that product does not factor into a face density times a
   * face volumetric flux.
   */
  virtual bool tracksOwnVolumetricFaceFlux() const { return false; }

  /// Whether the predictor mass flux interpolates the density-velocity product on internal faces
  const bool _momentum_interpolated_predictor_mass_flux;

  /// Interpolation method used for the pressure diffusion coefficient on faces
  const Moose::FV::InterpMethod _pressure_diffusion_interp_method;

private:
  /// The subset of the FaceInfo objects that actually cover the subdomains which the
  /// flow field is defined on. Cached for performance optimization.
  std::vector<const FaceInfo *> _flow_face_info;
};

template <typename VarType>
void
RhieChowMassFlux::checkBlocks(const VarType & var) const
{
  const auto & var_blocks = var.blockIDs();
  const auto & uo_blocks = blockIDs();

  // Error if this UO has any blocks that the variable does not
  std::set<SubdomainID> uo_blocks_minus_var_blocks;
  std::set_difference(uo_blocks.begin(),
                      uo_blocks.end(),
                      var_blocks.begin(),
                      var_blocks.end(),
                      std::inserter(uo_blocks_minus_var_blocks, uo_blocks_minus_var_blocks.end()));
  if (uo_blocks_minus_var_blocks.size() > 0)
    mooseError("Block restriction of interpolator user object '",
               this->name(),
               "' (",
               Moose::stringify(blocks()),
               ") includes blocks not in the block restriction of variable '",
               var.name(),
               "' (",
               Moose::stringify(var.blocks()),
               ")");

  // Get the blocks in the variable but not this UO
  std::set<SubdomainID> var_blocks_minus_uo_blocks;
  std::set_difference(var_blocks.begin(),
                      var_blocks.end(),
                      uo_blocks.begin(),
                      uo_blocks.end(),
                      std::inserter(var_blocks_minus_uo_blocks, var_blocks_minus_uo_blocks.end()));

  // For each block in the variable but not this UO, error if there is connection
  // to any blocks on the UO.
  for (auto & block_id : var_blocks_minus_uo_blocks)
  {
    const auto connected_blocks = _moose_mesh.getBlockConnectedBlocks(block_id);
    std::set<SubdomainID> connected_blocks_on_uo;
    std::set_intersection(connected_blocks.begin(),
                          connected_blocks.end(),
                          uo_blocks.begin(),
                          uo_blocks.end(),
                          std::inserter(connected_blocks_on_uo, connected_blocks_on_uo.end()));
    if (connected_blocks_on_uo.size() > 0)
      mooseError("Block restriction of interpolator user object '",
                 this->name(),
                 "' (",
                 Moose::stringify(uo_blocks),
                 ") doesn't match the block restriction of variable '",
                 var.name(),
                 "' (",
                 Moose::stringify(var_blocks),
                 ")");
  }
}
