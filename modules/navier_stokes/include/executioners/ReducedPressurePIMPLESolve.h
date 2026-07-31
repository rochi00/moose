//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PIMPLESolve.h"
#include "RhieChowMassFlux.h"

#include "libmesh/point.h"

#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class ConservativeSharpInterfaceRhieChowMassFlux;
class ConservativeSharpInterfaceVOFMULESCorrector;
class MassConservedLevelSetCorrector;
class LowMachDivergenceSource;
class LowMachEnthalpyNewtonState;
class LowMachImplicitMassFlux;

/**
 * PIMPLE solve object with explicit hooks for additional reduced-pressure /
 * sharp-interface face-flux predictors.
 */
class ReducedPressurePIMPLESolve : public PIMPLESolve
{
public:
  static InputParameters validParams();

  ReducedPressurePIMPLESolve(Executioner & ex);

  void commitAcceptedTimestepTransportHistory() const;

protected:
  void initialSetup() override;
  void preSolveSetup(const SolverParams & solver_params) override;
  void addIterationResiduals(ResidualStorage & residual_storage) override;
  void initializeSolveLoop(const SolverParams & solver_params) override;
  void preMomentumPressureIteration(ResidualStorage & residual_storage,
                                    const SolverParams & solver_params) override;
  bool shouldAssembleMomentumPredictorWithoutSolve() const override;
  bool shouldSolveEnergyAfterFlowLoop() const override;
  bool shouldSolveActiveScalarsAfterFlowLoop() const override;
  void finalizeSolve(const bool converged) override;
  void addMomentumPredictorExplicitForcing(const unsigned int system_i,
                                           NumericVector<Number> & rhs) override;
  bool shouldCopyMomentumNonlinearSolutionHistory() const override { return false; }
  void postPreparePressureCorrectorState(const bool subtract_updated_pressure) override;
  void postPublishPressureCorrectedState() override;

private:
  struct LowMachEnthalpyResidual
  {
    Real absolute;
    Real normalization;
    Real algebraic_scale;

    Real normalized() const { return absolute / normalization; }
    bool converged(const Real relative_tolerance) const
    {
      constexpr Real roundoff_multiplier = 100.0;
      return normalized() <= relative_tolerance ||
             absolute <=
                 roundoff_multiplier * std::numeric_limits<Real>::epsilon() * algebraic_scale;
    }
  };

  bool startupPressureInitializationEnabled() const;
  bool shouldRunStartupInitialization() const;
  bool solvesVolumeFraction() const;
  void initializeStartupPressureField(const SolverParams & solver_params);
  void resetVOFTransportStateForNewSolve() const;
  void initializeConsistentStartupState();
  void captureLowMachOuterState();
  Real updateLowMachOuterState();
  void applyLevelSetAitkenRelaxation(unsigned int system_index,
                                     LinearSystem & system,
                                     MassConservedLevelSetCorrector & corrector);
  void commitAcceptedVOFTransportHistoryIfNeeded() const;
  void advanceOuterIterationHistories();
  void solveImplicitLowMachMassEquation(ResidualStorage & residual_storage,
                                        const SolverParams & solver_params);
  void solveImplicitLowMachEnthalpyEquation(ResidualStorage & residual_storage,
                                            const SolverParams & solver_params);
  LowMachEnthalpyResidual assembleLowMachEnthalpyResidual();
  void validateLowMachEnthalpyOperators() const;
  void solveVolumeFractionBeforeFlowCorrection(ResidualStorage & residual_storage,
                                               const SolverParams & solver_params);
  void prepareVOFTransportStateForOuterIteration() const;
  void adoptPublishedVOFTransportState() const;
  void
  storeActiveScalarResiduals(ResidualStorage & residual_storage,
                             const std::vector<std::pair<unsigned int, Real>> & vf_residuals) const;
  unsigned int computeVolumeFractionSubcycles() const;
  void synchronizeSystemState(LinearSystem & system) const;
  void setPreviousNewtonToCurrent(LinearSystem & system) const;
  void advanceVolumeFractionSubcycleOldState(LinearSystem & system) const;
  void setProblemSubcycleTime(const unsigned int subcycle,
                              const Real subcycle_dt,
                              const Real global_time_old);
  std::vector<std::pair<unsigned int, Real>> solveVolumeFractionSystems();
  std::pair<unsigned int, Real> solveOneVolumeFractionSystem(const unsigned int i,
                                                             const unsigned int num_subcycles,
                                                             const Real subcycle_dt,
                                                             const Real global_dt,
                                                             const Real global_time_old);
  std::pair<unsigned int, Real>
  runOneVolumeFractionSubcycle(const unsigned int i,
                               LinearSystem & system,
                               ConservativeSharpInterfaceVOFMULESCorrector & corrector,
                               const unsigned int subcycle,
                               const Real subcycle_dt,
                               const Real global_dt,
                               const Real global_time_old);
  void finalizeVolumeFractionTransportState();
  ConservativeSharpInterfaceVOFMULESCorrector *
  sharpInterfaceVOFCorrector(const SolverSystemName & system_name) const;
  MassConservedLevelSetCorrector *
  massConservedLevelSetCorrector(const SolverSystemName & system_name) const;

  std::pair<unsigned int, Real> correctStartupContinuityOnce(const SolverParams & solver_params);

  ConservativeSharpInterfaceRhieChowMassFlux * sharpInterfaceRC() const;

  const unsigned int _volume_fraction_subcycles;
  const Real _volume_fraction_max_courant;
  std::string _startup_pressure_initialization;
  const unsigned int _startup_flux_corrections;

  const SolverSystemName _low_mach_mass_system_name;
  const bool _has_low_mach_mass_system;
  const unsigned int _low_mach_mass_system_number;
  LinearSystem * const _low_mach_mass_system;
  LowMachImplicitMassFlux * _low_mach_mass_flux;
  Moose::PetscSupport::PetscOptions _low_mach_mass_petsc_options;
  SIMPLESolverConfiguration _low_mach_mass_linear_control;
  const Real _low_mach_mass_l_abs_tol;
  const Real _low_mach_mass_absolute_tolerance;
  std::size_t _low_mach_mass_residual_index = Moose::invalid_size_t;
  const unsigned int _low_mach_minimum_outer_iterations;
  const Real _low_mach_fixed_point_tolerance;
  const Real _low_mach_face_flux_absolute_tolerance;
  const Real _low_mach_level_set_aitken_min_relaxation;
  const Real _low_mach_level_set_aitken_max_relaxation;
  const Real _low_mach_level_set_aitken_max_growth;
  const Real _low_mach_level_set_aitken_projection_amplification;
  std::size_t _low_mach_fixed_point_residual_index = Moose::invalid_size_t;
  std::vector<std::unique_ptr<NumericVector<Number>>> _low_mach_previous_outer_solutions;
  std::unordered_map<dof_id_type, Real> _low_mach_previous_outer_face_flux;
  std::vector<std::unique_ptr<NumericVector<Number>>> _level_set_previous_aitken_residuals;
  std::vector<Real> _level_set_aitken_relaxation;
  std::vector<Real> _level_set_raw_fixed_point_update;

  LowMachEnthalpyNewtonState * _low_mach_enthalpy_state;
  LowMachDivergenceSource * _low_mach_divergence_source;
  const unsigned int _low_mach_enthalpy_max_iterations;
  const Real _low_mach_enthalpy_tolerance;
  const Real _low_mach_enthalpy_residual_tolerance;
  const Real _low_mach_enthalpy_l_tol;
  std::size_t _low_mach_liquid_fraction_residual_index = Moose::invalid_size_t;
};
