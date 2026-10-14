//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "RhieChowMassFlux.h"
#include "DerivedFaceFunctor.h"
#include "PetscVectorReader.h"

class ElemInfo;
class LinearSystem;

/**
 * Rhie-Chow face-flux provider extended with additional reduced-pressure source
 * fluxes for large-density-ratio diffuse-interface coupling.
 *
 * Sign convention:
 * - The face functors published here are PRESSURE-EQUATION SOURCE FLUXES in the
 *   same sign convention as the stock HbyA functor used by LinearFVDivergence.
 * - The physical face mass-flux correction is therefore the negative of the
 *   published source functor. This mirrors the stock MOOSE relation
 *     face_mass_flux = -HbyA + pressure_diffusion_flux.
 */
class ConservativeDiffuseInterfaceRhieChowMassFluxBase : public RhieChowMassFlux
{
public:
  struct DiffuseFaceOperatorState
  {
    Real face_rho = 0.0;
    Real normal_raw_ainv = 0.0;
    RealVectorValue face_normal;
    RealVectorValue face_raw_ainv;
  };

  static InputParameters validParams();

  ConservativeDiffuseInterfaceRhieChowMassFluxBase(const InputParameters & params);

  void meshChanged() override;
  void initialSetup() override;
  void initialize() override;
  void beginMomentumPredictorAssembly(const unsigned int system_i) override;
  void momentumPredictorAssembled(const unsigned int system_i,
                                  libMesh::SparseMatrix<Number> & matrix,
                                  const NumericVector<Number> & rhs,
                                  const Real relaxation_parameter) override;
  /// The predictor is assembled without the pressure and body forces, which
  /// addMomentumPredictorExplicitForcing supplies from the reduced pressure, gravity and
  /// capillarity; H/A then comes from the operator as assembled, treating the boundary
  /// coefficients as OpenFOAM's fvMatrix::A() and H() do.
  bool splitMomentumPredictorOperator() const override { return true; }
  void computeHbyA(const bool verbose) override;
  void cachePressureEquationFlux() override;
  void computeFaceMassFlux() override;
  void computeCellVelocity() override;
  Real getMassFlux(const FaceInfo & fi) const override;
  Real getVolumetricFaceFlux(const FaceInfo & fi) const override;
  void commitAcceptedTimestepTransportHistory();
  void freezeLaggedCorrectedFaceFlux(const bool use_previous_timestep_flux = false);

  /// Update the additional pressure-equation source-flux functors before the pressure solve.
  void updateAdditionalPressureFluxFunctors(const bool verbose);
  /// Update boundary pressure gradients from the current predictor and boundary velocity state.
  void updatePressureBoundaryNormalGradients();
  /// Refresh boundary-face velocity values from the active FV velocity BC objects.
  void updateVelocityBoundaryState();
  /// Maximum face-flux Courant number over active flow cells for the supplied timestep.
  Real maxCourant(const Real dt) const;

  void computePressureCorrectedCellVelocity();

  /**
   * While the startup continuity projection runs, the pressure equation carries no transient,
   * hydrostatic or capillary source and the predictor force carries no hydrostatic part: that
   * projection repairs the face flux alone, and leaves the user's reduced pressure untouched.
   */
  void setStartupProjectionActive(const bool active) { _startup_projection_active = active; }
  bool startupProjectionActive() const { return _startup_projection_active; }
  RealVectorValue
  reducedPressureMomentumPredictorForceDensity(const ElemInfo & elem_info,
                                               const Moose::StateArg & time_arg) const;
  Real predictorVelocityComponent(const ElemInfo & elem_info, const unsigned int component) const;
  Real vofRhoPhiIntegrated(const FaceInfo & fi) const;

protected:
  using FaceScalarField = FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>>;
  struct MomentumPredictorOperator;

  void computePredictorOperator(const unsigned int system_i,
                                NumericVector<Number> & HbyA,
                                NumericVector<Number> & Ainv) override;
  /// Finalize all component operators as one vector momentum-predictor object.
  void finalizeCachedMomentumPredictorOperators();
  /// Access the momentum-predictor operator for a component.
  const MomentumPredictorOperator *
  cachedMomentumPredictorOperator(const unsigned int system_i) const;

  /// Target boundary flux of the pressure equation: the prescribed boundary volumetric flux
  Real pressureBoundaryTargetFlux(const FaceInfo * fi, const Moose::StateArg & time_arg) const;
  /// Normal diffusion coefficient of the pressure equation on a boundary face
  Real pressureBoundaryNormalAinv(const FaceInfo * fi) const;
  bool pressureCorrectionFluxIsIntegrated() const override { return true; }
  Real boundaryVelocityValue(const FaceInfo * fi,
                             const unsigned int component,
                             const Moose::StateArg & time_arg) const override;
  /// Compute the target volumetric boundary flux U_b.n at a face.
  Real boundaryVolumetricFluxTarget(const FaceInfo * fi, const Moose::StateArg & time_arg) const;
  /// Compute the normal component of a face diagonal Ainv vector.
  Real computeFaceNormalAinv(const RealVectorValue & face_ainv,
                             const RealVectorValue & face_normal) const;
  /// Interpolate per-component momentum-system cell readers to a face vector.
  RealVectorValue
  interpolateMomentumVectorReadersToFace(const FaceInfo * fi,
                                         const std::vector<PetscVectorReader> & readers) const;
  Real pressureFaceNormalAinv(const FaceInfo * fi, const Moose::StateArg & time_arg) const;

  void initializeAdditionalPressureFluxStorage();
  void writePressureCorrectedVelocityToMomentumSolution(const Moose::StateArg & time_arg);
  Real correctedTransportPhi(const FaceInfo & fi) const;
  void cacheCurrentCorrectedVolumetricFlux(const Real degenerate_normal_pressure_ainv_tol = 0.0);
  /// The face-normal velocity correction the pressure equation produced, divided by the face
  /// coefficient it was produced through, for the cell-centered velocity writeback.
  Real writebackReconstructionScalar(const FaceInfo * fi,
                                     const Moose::StateArg & time_arg,
                                     const Real degenerate_normal_pressure_ainv_tol = 0.0) const;
  Real maxPressureFaceNormalAinv(const Moose::StateArg & time_arg) const;
  Real transportMassFluxDensityFromVolumetricPhi(const FaceInfo * fi,
                                                 const Real volumetric_phi,
                                                 const Moose::StateArg & time_arg) const;

  Real interpolateFaceDensity(const FaceInfo * fi, const Moose::StateArg & time_arg) const;
  Real predictorFaceDensity(const FaceInfo * fi, const Moose::StateArg & time_arg) const;
  Real facePhysicalVelocityComponent(const FaceInfo * fi,
                                     const unsigned int component,
                                     const Moose::StateArg & time_arg) const;
  Real interpolatedPhysicalFaceFlux(const FaceInfo * fi, const Moose::StateArg & time_arg) const;

  DiffuseFaceOperatorState
  buildDiffuseFaceOperatorState(const FaceInfo * fi,
                                const Moose::StateArg & time_arg,
                                const std::vector<PetscVectorReader> & raw_ainv_readers) const;

  Real evaluateFaceScalarFunctor(const Moose::Functor<Real> * functor,
                                 const FaceInfo * fi,
                                 const Moose::StateArg & time_arg) const;

  Real computeDefaultTransientProjectionVolumetricFlux(const FaceInfo * fi,
                                                       const Moose::StateArg & time_arg,
                                                       const DiffuseFaceOperatorState & state);
  Real computeFaceNormalDensityGradient(const FaceInfo * fi,
                                        const Moose::StateArg & time_arg) const;
  Real computeFaceNormalPressureGradient(const FaceInfo * fi,
                                         const Moose::StateArg & time_arg) const;
  bool useExplicitHydrostaticPredictorForce() const;
  Real
  computeDiscreteHydrostaticPredictorFaceNormalForceDensity(const FaceInfo * fi,
                                                            const Moose::StateArg & time_arg) const;
  bool surfaceTensionEnabled() const;
  /// The published capillary force, zeroed while the startup projection suppresses flux sources.
  RealVectorValue surfaceTensionForceDensity(const ElemInfo & elem_info,
                                             const Moose::StateArg & time_arg) const;
  Real surfaceTensionFaceNormalForceDensity(const FaceInfo * fi,
                                            const Moose::StateArg & time_arg) const;
  void updatePressureCoupledVelocityCorrectionFaceField(const Moose::StateArg & time_arg);
  RealVectorValue
  reconstructFaceNormalScalarToCellVector(const ElemInfo * elem_info,
                                          const Moose::StateArg & time_arg,
                                          const FaceScalarField & scalar_field) const;
  RealVectorValue
  reconstructBasePressureCoupledCellVelocityDelta(const ElemInfo * elem_info,
                                                  const Moose::StateArg & time_arg) const;
  virtual RealVectorValue
  reconstructPressureCoupledCellVelocityDelta(const ElemInfo * elem_info,
                                              const Moose::StateArg & time_arg) const;
  bool useConstrainedBoundaryPredictorState(const FaceInfo * fi) const override;

  /**
   * Explicit predictor-side face flux entering the pressure equation in the native
   * pressure-correction sign convention. The accepted transport-face flux is therefore formed as
   * -_phiHbyA_flux + _pressure_equation_flux.
   */
  DerivedFaceFunctor<Real> _phiHbyA_flux;
  /// Explicit face-force contribution entering the pressure corrector.
  FaceScalarField _phig_flux;
  /// Face diffusion tensor of the pressure equation: the raw face rAUf.
  FaceCenteredMapFunctor<RealVectorValue, std::unordered_map<dof_id_type, RealVectorValue>>
      _pressure_Ainv;
  /// Boundary-normal pressure gradient the fixed-flux pressure boundaries impose
  FaceScalarField _pressure_boundary_normal_gradient;
  /// Predictor-side face flux before the explicit face force is removed, stored in the internal
  /// pressure-correction convention. The accepted transport flux uses the negated value.
  std::unordered_map<dof_id_type, Real> _pressure_predictor_base_flux;
  /// Boundary-face velocity values refreshed from the current FV velocity BCs.
  std::vector<std::unordered_map<dof_id_type, Real>> _boundary_velocity_face_values;
  /// Whether _boundary_velocity_face_values matches the latest cell-centered velocity.
  bool _velocity_boundary_state_valid = false;
  /// Whether _pressure_boundary_normal_gradient matches the latest predictor state.
  bool _pressure_boundary_normal_gradient_valid = false;
  /// Whether _phiHbyA_flux/_phig_flux match the latest predictor state.
  bool _pressure_predictor_face_state_valid = false;

  DerivedFaceFunctor<Real> _corrected_face_phi;
  FaceScalarField _previous_timestep_corrected_face_phi;
  FaceScalarField _lagged_corrected_face_phi;
  FaceScalarField _pressure_coupled_cell_reconstruction_scalar;

  const RealVectorValue _gravity;
  const Point _reference_pressure_point;

  const MooseFunctorName _vof_rho_phi_name;
  const MooseFunctorName _surface_tension_force_name;

  const Moose::Functor<Real> * _vof_rho_phi;
  const Moose::Functor<RealVectorValue> * _surface_tension_force = nullptr;
  bool _startup_projection_active = false;
  bool _pressure_coupled_velocity_correction_valid = false;

  /// Predictor operators for each momentum component.
  std::vector<std::unique_ptr<MomentumPredictorOperator>> _momentum_predictor_operators;

  /// Whether every component operator has been assembled and finalized since the last begin
  bool _predictor_ready = false;
};
