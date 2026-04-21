//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// Moose includes
#include "SolveObject.h"
#include "UserObjectInterface.h"
#include "PetscSupport.h"
#include "SolverParams.h"
#include "SegregatedSolverUtils.h"
#include "RhieChowMassFluxMultiPhase.h"

class LinearFVMultiPhaseFractionAdvection;

// Libmesh includes
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/petsc_matrix.h"
#include "libmesh/equation_systems.h"

/**
 * Solver configuration class used with the linear solvers in a SIMPLE solver.
 */
class SolverConfigurationMultiPhase : public libMesh::SolverConfiguration
{
  /**
   * Override this to make sure the PETSc options are not overwritten in the linear solver
   */
  virtual void configure_solver() override {}
};

/**
 * Solve class serving as a base class for the two SIMPLE solvers that operate with
 * different assembly algorithms. Includes base routines and variables for the coupling of
 * momentum and pressure.
 */
class SolveBaseMultiPhase : public SolveObject, public UserObjectInterface
{
public:
  SolveBaseMultiPhase(Executioner & ex);

  static InputParameters validParams();

  virtual void setInnerSolve(SolveObject &) override
  {
    mooseError("Cannot set inner solve object for solves that inherit from SolveBaseMultiPhase");
  }

  /// Fetch the Rhie Chow user object that is reponsible for determining face
  /// velocities and mass flux
  virtual void linkRhieChowUserObjects();

  /// Setup pressure pin if there is need for one
  void setupPressurePin();

  /// Check if the user defined time kernels
  virtual void checkIntegrity() {}

  /**
   * Performs the momentum pressure coupling.
   * @return True if solver is converged.
   */
  virtual bool solve() override;

  /// Return pointers to the systems which are solved for within this object
  const std::vector<LinearSystem *> systemsToSolve() const { return _systems_to_solve; }

protected:
  void checkDependentParameterError(const std::string & main_parameter,
                                    const std::vector<std::string> & dependent_parameters,
                                    const bool should_be_defined);

  /// Number of phases in the solve
  const unsigned int _number_of_phases;

  /// Solve a momentum predictor step with a fixed pressure field
  /// @return A vector of (number of linear iterations, normalized residual norm) pairs for
  /// the momentum equations. The length of the vector equals the dimensionality of
  /// the domain.
  std::vector<std::pair<unsigned int, Real>> solveMomentumPredictor(const unsigned int phase_number);

  /// Solve a pressure corrector step.
  /// @return The number of linear iterations and the normalized residual norm of
  /// the pressure equation.
  std::pair<unsigned int, Real> solvePressureCorrector();

  /// Computes new velocity field based on computed pressure gradients
  /// @param phase_number The number of the phase for which the problem is being solved
  /// @param subtract_updated_pressure If we need to subtract the updated
  /// pressure gradient from the right hand side of the system
  /// @param recompute_face_mass_flux If we want to recompute the face flux too
  /// @param solver_params Dummy solver parameter object for the linear solve
  virtual std::pair<unsigned int, Real> correctVelocities(const bool subtract_updated_pressure,
                                                          const bool recompute_face_mass_flux,
                                                          const SolverParams & solver_params);

  /// Solve an equation which contains an advection term that depends
  /// on the solution of the segregated Navier-Stokes equations.
  /// @param system_num The number of the system which is solved
  /// @param system Reference to the system which is solved
  /// @param relaxation_factor The relaxation factor for matrix relaxation
  /// @param solver_config The solver configuration object for the linear solve
  /// @param abs_tol The scaled absolute tolerance for the linear solve
  /// @param relax_fields (optional) A boolean flag to indicate whether to relax fields during the solve. Default value is false.
  /// @param field_relaxation (optional) The relaxation factor for fields if relax_fields is true. Default value is 1.0.
  /// @return The normalized residual norm of the equation.
  std::pair<unsigned int, Real> solveAdvectedSystem(const unsigned int system_num,
                                                    LinearSystem & system,
                                                    const Real relaxation_factor,
                                                    libMesh::SolverConfiguration & solver_config,
                                                    const Real abs_tol,
                                                    const bool relax_fields = false,
                                                    const Real field_relaxation = 1.0);

  // ************************ Momentum Eq Variables ************************ //

  /// The names of the momentum systems.
  const std::vector<std::vector<std::string>> & _momentum_system_names;

  /// Options for the linear solver of the momentum equation
  SolverConfigurationMultiPhase _momentum_linear_control;

  /// Absolute linear tolerance for the momentum equation(s). We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _momentum_l_abs_tol;

  /// Options which hold the petsc settings for the momentum equation
  Moose::PetscSupport::PetscOptions _momentum_petsc_options;

  /// The user-defined relaxation parameter for the momentum equation
  const Real _momentum_equation_relaxation;

  // ************************ Pressure Eq Variables ************************ //

  /// The name of the pressure system
  const SolverSystemName & _pressure_system_name;

  /// Options for the linear solver of the pressure equation
  SolverConfigurationMultiPhase _pressure_linear_control;

  /// Absolute linear tolerance for the pressure equation. We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _pressure_l_abs_tol;

  /// Options which hold the petsc settings for the pressure equation
  Moose::PetscSupport::PetscOptions _pressure_petsc_options;

  /// The user-defined relaxation parameter for the pressure variable
  const Real _pressure_variable_relaxation;

  /// If the pressure needs to be pinned
  const bool _pin_pressure;

  /// The value we want to enforce for pressure
  const Real _pressure_pin_value;

  /// The dof ID where the pressure needs to be pinned
  dof_id_type _pressure_pin_dof;

  // ************************ Energy Eq Variables ************************** //

  /// Boolean for easy check if a fluid energy system shall be solved or not
  const bool _has_energy_system;

  /// The names of the fluid energy scalar systems
  const std::vector<SolverSystemName> * _energy_system_names;

  /// The user-defined relaxation parameter for the energy equation
  const Real _energy_equation_relaxation;

  /// Options which hold the petsc settings for the fluid energy equation
  Moose::PetscSupport::PetscOptions _energy_petsc_options;

  /// Options for the linear solver of the energy equation
  SolverConfigurationMultiPhase _energy_linear_control;

  /// Absolute linear tolerance for the energy equations. We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _energy_l_abs_tol;

  // ************************ Solid Energy Eq Variables *********************** //

  /// Boolean for easy check if a solid energy system shall be solved or not
  const bool _has_solid_energy_system;

  /// The name of the solid energy scalar systems
  const SolverSystemName _solid_energy_system_names;

  /// The user-defined relaxation parameter for the solid energy equation
  const Real _solid_energy_equation_relaxation;

  /// Options which hold the petsc settings for the fluid energy equation
  Moose::PetscSupport::PetscOptions _solid_energy_petsc_options;

  /// Options for the linear solver of the energy equation
  SolverConfigurationMultiPhase _solid_energy_linear_control;

  /// Absolute linear tolerance for the energy equations. We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _solid_energy_l_abs_tol;

  // ************************ Phase Transport System ************************** //

  /// The names of the phases
  const std::vector<SolverSystemName> & _phase_system_names;

  /// Number of phases solved for the void fraction
  const unsigned int _number_of_solved_phases;

  /// The user-defined relaxation parameter for the phase equation
  const Real _phase_equation_relaxation;

  /// Number of internal MULES iterations to perform
  const unsigned int _MULES_iterations;

  /// Number of alpha sub-cycles per time step
  const unsigned int _n_alpha_subcycles;

  /// Options which hold the petsc settings for the phase equation
  Moose::PetscSupport::PetscOptions _phase_petsc_options;

  /// Options for the linear solver of the phase equation
  SolverConfigurationMultiPhase _phase_linear_control;

  /// Absolute linear tolerance for the phase equations. We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _phase_l_abs_tol;

  // ************************ Passive Scalar Variables ************************ //

  /// The names of the passive scalar systems
  const std::vector<SolverSystemName> & _passive_scalar_system_names;

  /// Boolean for easy check if a passive scalar systems shall be solved or not
  const bool _has_passive_scalar_systems;

  // The number(s) of the system(s) corresponding to the passive scalar equation(s)
  std::vector<unsigned int> _passive_scalar_system_numbers;

  /// The user-defined relaxation parameter(s) for the passive scalar equation(s)
  const std::vector<Real> _passive_scalar_equation_relaxation;

  /// Options which hold the petsc settings for the passive scalar equation(s)
  Moose::PetscSupport::PetscOptions _passive_scalar_petsc_options;

  /// Options for the linear solver of the passive scalar equation(s)
  SolverConfigurationMultiPhase _passive_scalar_linear_control;

  /// Absolute linear tolerance for the passive scalar equation(s). We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _passive_scalar_l_abs_tol;

  // ************************ Turbulence Variables ************************ //

  /// The names of the turbulence systems
  const std::vector<std::vector<std::string>> * _turbulence_system_names;

  /// Boolean for easy check if a turbulence scalar systems shall be solved or not
  const bool _has_turbulence_systems;

  // The number(s) of the system(s) corresponding to the turbulence equation(s)
  std::vector<std::vector<unsigned int>> _turbulence_system_numbers;

  /// The user-defined relaxation parameter(s) for the turbulence equation(s)
  const std::vector<Real> _turbulence_equation_relaxation;

  /// The user-defined relaxation parameter(s) for the turbulence field(s)
  std::vector<Real> _turbulence_field_relaxation;

  /// The user-defined lower limit for turbulent quantities e.g. k, eps/omega, etc..
  std::vector<Real> _turbulence_field_min_limit;

  /// Options which hold the petsc settings for the turbulence equation(s)
  Moose::PetscSupport::PetscOptions _turbulence_petsc_options;

  /// Options for the linear solver of the turbulence equation(s)
  SolverConfigurationMultiPhase _turbulence_linear_control;

  /// Absolute linear tolerance for the turbulence equation(s). We need to store this, because
  /// it needs to be scaled with a representative flux.
  const Real _turbulence_l_abs_tol;

  // ************************ Iteration control **************************** //

  /// The user-defined absolute tolerance for determining the convergence in momentum
  const Real _momentum_absolute_tolerance;

  /// The user-defined absolute tolerance for determining the convergence in pressure
  const Real _pressure_absolute_tolerance;

  /// The user-defined absolute tolerance for determining the convergence in energy
  const Real _energy_absolute_tolerance;

  /// The user-defined absolute tolerance for determining the convergence in solid energy
  const Real _solid_energy_absolute_tolerance;

  /// The user-defined absolute tolerance for determining the convergence in the phase solutions
  const Real _phase_absolute_tolerance;

  /// The user-defined absolute tolerance for determining the convergence in passive scalars
  const std::vector<Real> _passive_scalar_absolute_tolerance;

  /// The user-defined absolute tolerance for determining the convergence turbulence variables
  const std::vector<Real> _turbulence_absolute_tolerance;

  /// The maximum number of momentum-pressure iterations
  const unsigned int _num_iterations;

  /// If solve should continue if maximum number of iterations is hit
  const bool _continue_on_max_its;

  // ************************ Other Variables ****************************** //

  /// Debug parameter which allows printing the coupling and solution vectors/matrices
  const bool _print_fields;

  // ********************** System Storage and Numbering ****************** //

  /// Solve an equation which contains the solid energy conservation.
  std::pair<unsigned int, Real> solveSolidEnergy();

  /// The number(s) of the system(s) corresponding to the momentum equation(s)
  std::vector<std::vector<unsigned int>> _momentum_system_numbers;

  /// Pointer(s) to the system(s) corresponding to the momentum equation(s)
  std::vector<std::vector<LinearSystem *>> _momentum_systems;

  /// The number of the system corresponding to the pressure equation
  const unsigned int _pressure_sys_number;

  /// Reference to the nonlinear system corresponding to the pressure equation
  LinearSystem & _pressure_system;

  /// The number of the system corresponding to the energy equation
  std::vector<unsigned int> _energy_system_numbers;

  /// Pointer to the nonlinear system corresponding to the fluid energy equation
  std::vector<LinearSystem *> _energy_systems;

  /// The number of the system corresponding to the solid energy equation
  const unsigned int _solid_energy_sys_number;

  /// Pointer to the nonlinear system corresponding to the solid energy equation
  LinearSystem * _solid_energy_system;

  /// The number of the system corresponding to the energy equation
  std::vector<unsigned int> _phase_system_numbers;

  /// Pointer to the nonlinear system corresponding to the fluid energy equation
  std::vector<LinearSystem *> _phase_systems;

  /// Pointer(s) to the system(s) corresponding to the passive scalar equation(s)
  std::vector<LinearSystem *> _passive_scalar_systems;

  /// Pointer(s) to the system(s) corresponding to the turbulence equation(s)
  std::vector<std::vector<LinearSystem *>> _turbulence_systems;

  /// Pointer to the segregated RhieChow interpolation object
  std::vector<RhieChowMassFluxMultiPhase *> _rc_uo;

  /// Shortcut to every linear system that we solve for here
  std::vector<LinearSystem *> _systems_to_solve;

  /// Interface sharpening variables
  const bool _enforce_phase_sum;
  const bool _activate_interface_shapening;
  const MooseEnum _shapening_type;
  const Real _smoothing_constant;

  /// Pointers to the alpha advection kernels (for semi-implicit MULES correction)
  std::vector<LinearFVMultiPhaseFractionAdvection *> _alpha_advection_kernels;

};
