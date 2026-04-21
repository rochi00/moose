//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolveBaseMultiPhase.h"
#include "FEProblem.h"
#include "SegregatedSolverUtils.h"
#include "LinearSystem.h"
#include "LinearFVMultiPhaseFractionAdvection.h"
#include "LinearFVFluxKernel.h"
#include "Attributes.h"

using namespace libMesh;

InputParameters
SolveBaseMultiPhase::validParams()
{
  InputParameters params = emptyInputParameters();
  params.addRequiredParam<std::vector<UserObjectName>>("rhie_chow_user_objects", "The rhie-chow user-objects");

  /*
   * Define the number of phases
   */
  params.addRangeCheckedParam<unsigned int>(
      "number_of_phases",
      1,
      "0 < number_of_phases",
      "The number of phases in the problem.");

  /*
   * The names of the different systems in the segregated solver
   */
  params.addRequiredParam<std::vector<std::vector<std::string>>>(
      "momentum_systems", "The solver system(s) for the momentum equation(s) for all phases.");
  params.addRequiredParam<SolverSystemName>("pressure_system",
                                            "The solver system for the pressure equation.");
  params.addParam<std::vector<SolverSystemName>>("energy_systems", "The solver system for the energy equation.");
  params.addParam<SolverSystemName>("solid_energy_system",
                                    "The solver system for the solid energy equation.");
  params.addParam<std::vector<SolverSystemName>>("phase_systems", "The solver system for the phase equations.");
  params.addParam<std::vector<SolverSystemName>>(
      "passive_scalar_systems", {}, "The solver system for each scalar advection equation.");
  params.addParam<std::vector<std::vector<std::string>>>(
      "turbulence_systems", "The solver system for each surrogate turbulence equation.");

  /*
   * Parameters to control the solution of the momentum equation
   */

  params.addRangeCheckedParam<Real>(
      "momentum_equation_relaxation",
      0.7,
      "0.0<momentum_equation_relaxation<=1.0",
      "The relaxation which should be used for the momentum equation. (=1 for no relaxation, "
      "diagonal dominance will still be enforced)");

  params.addParam<MultiMooseEnum>("momentum_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the momentum equation");
  params.addParam<MultiMooseEnum>("momentum_petsc_options_iname",
                                  Moose::PetscSupport::getCommonPetscKeys(),
                                  "Names of PETSc name/value pairs for the momentum equation");
  params.addParam<std::vector<std::string>>(
      "momentum_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "momentum equation");

  params.addRangeCheckedParam<Real>(
      "momentum_absolute_tolerance",
      1e-5,
      "0.0<momentum_absolute_tolerance",
      "The absolute tolerance on the normalized residual of the momentum equation.");

  params.addRangeCheckedParam<Real>("momentum_l_tol",
                                    1e-5,
                                    "0.0<=momentum_l_tol & momentum_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the momentum equation.");
  params.addRangeCheckedParam<Real>("momentum_l_abs_tol",
                                    1e-50,
                                    "0.0<momentum_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the momentum equation.");
  params.addParam<unsigned int>(
      "momentum_l_max_its",
      10000,
      "The maximum allowed iterations in the linear solver of the momentum equation.");

  params.addParamNamesToGroup(
      "momentum_equation_relaxation momentum_petsc_options momentum_petsc_options_iname "
      "momentum_petsc_options_value momentum_petsc_options_value momentum_absolute_tolerance "
      "momentum_l_tol momentum_l_abs_tol momentum_l_max_its momentum_systems",
      "Momentum Equation");

  /*
   * Parameters to control the solution of the pressure equation
   */
  params.addRangeCheckedParam<Real>(
      "pressure_variable_relaxation",
      0.3,
      "0.0<pressure_variable_relaxation<=1.0",
      "The relaxation which should be used for the pressure variable (=1 for no relaxation).");

  params.addParam<MultiMooseEnum>("pressure_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the pressure equation");
  params.addParam<MultiMooseEnum>("pressure_petsc_options_iname",
                                  Moose::PetscSupport::getCommonPetscKeys(),
                                  "Names of PETSc name/value pairs for the pressure equation");
  params.addParam<std::vector<std::string>>(
      "pressure_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "pressure equation");

  params.addRangeCheckedParam<Real>(
      "pressure_absolute_tolerance",
      1e-5,
      "0.0<pressure_absolute_tolerance",
      "The absolute tolerance on the normalized residual of the pressure equation.");

  params.addRangeCheckedParam<Real>("pressure_l_tol",
                                    1e-5,
                                    "0.0<=pressure_l_tol & pressure_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the pressure equation.");
  params.addRangeCheckedParam<Real>("pressure_l_abs_tol",
                                    1e-10,
                                    "0.0<pressure_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the pressure equation.");
  params.addParam<unsigned int>(
      "pressure_l_max_its",
      10000,
      "The maximum allowed iterations in the linear solver of the pressure equation.");

  params.addParamNamesToGroup(
      "pressure_variable_relaxation pressure_petsc_options pressure_petsc_options_iname "
      "pressure_petsc_options_value pressure_petsc_options_value pressure_absolute_tolerance "
      "pressure_l_tol pressure_l_abs_tol pressure_l_max_its pressure_system",
      "Pressure Equation");

  /*
   * Pressure pin parameters for enclosed flows
   */

  params.addParam<bool>(
      "pin_pressure", false, "If the pressure field needs to be pinned at a point.");
  params.addParam<Real>(
      "pressure_pin_value", 0.0, "The value which needs to be enforced for the pressure.");
  params.addParam<Point>("pressure_pin_point", "The point where the pressure needs to be pinned.");

  params.addParamNamesToGroup("pin_pressure pressure_pin_value pressure_pin_point", "Pressure Pin");

  params.addParam<bool>(
      "print_fields",
      false,
      "Use this to print the coupling and solution fields and matrices throughout the iteration.");

  /*
   * Parameters to control the solution of the energy equation
   */

  params.addRangeCheckedParam<Real>(
      "energy_equation_relaxation",
      0.9,
      "0.0<energy_equation_relaxation<=1.0",
      "The relaxation which should be used for the energy equation. (=1 for no relaxation, "
      "diagonal dominance will still be enforced)");

  params.addParam<MultiMooseEnum>("energy_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the energy equation");
  params.addParam<MultiMooseEnum>("energy_petsc_options_iname",
                                  Moose::PetscSupport::getCommonPetscKeys(),
                                  "Names of PETSc name/value pairs for the energy equation");
  params.addParam<std::vector<std::string>>(
      "energy_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "energy equation");

  params.addRangeCheckedParam<Real>(
      "energy_absolute_tolerance",
      1e-5,
      "0.0<energy_absolute_tolerance",
      "The absolute tolerance on the normalized residual of the energy equation.");

  params.addRangeCheckedParam<Real>("energy_l_tol",
                                    1e-5,
                                    "0.0<=energy_l_tol & energy_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the energy equation.");
  params.addRangeCheckedParam<Real>("energy_l_abs_tol",
                                    1e-10,
                                    "0.0<energy_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the energy equation.");
  params.addRangeCheckedParam<unsigned int>(
      "energy_l_max_its",
      10000,
      "0<energy_l_max_its",
      "The maximum allowed iterations in the linear solver of the energy equation.");

  params.addParamNamesToGroup(
      "energy_equation_relaxation energy_petsc_options energy_petsc_options_iname "
      "energy_petsc_options_value energy_petsc_options_value energy_absolute_tolerance "
      "energy_l_tol energy_l_abs_tol energy_l_max_its",
      "Energy Equation");

  /*
   * Parameters to control the solution of the solid energy equation
   */

  params.addRangeCheckedParam<Real>(
      "solid_energy_equation_relaxation",
      0.9,
      "0.0<solid_energy_equation_relaxation<=1.0",
      "The relaxation which should be used for the solid energy equation. (=1 for no relaxation, "
      "diagonal dominance will still be enforced)");

  params.addParam<MultiMooseEnum>("solid_energy_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the solid energy equation");
  params.addParam<MultiMooseEnum>("solid_energy_petsc_options_iname",
                                  Moose::PetscSupport::getCommonPetscKeys(),
                                  "Names of PETSc name/value pairs for the solid energy equation");
  params.addParam<std::vector<std::string>>(
      "solid_energy_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "solid energy equation");

  params.addRangeCheckedParam<Real>(
      "solid_energy_absolute_tolerance",
      1e-5,
      "0.0<solid_energy_absolute_tolerance",
      "The absolute tolerance on the normalized residual of the solid energy equation.");

  params.addRangeCheckedParam<Real>("solid_energy_l_tol",
                                    1e-5,
                                    "0.0<=solid_energy_l_tol & solid_energy_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the solid energy equation.");

  params.addRangeCheckedParam<Real>("solid_energy_l_abs_tol",
                                    1e-10,
                                    "0.0<solid_energy_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the solid energy equation.");
  params.addRangeCheckedParam<unsigned int>(
      "solid_energy_l_max_its",
      10000,
      "0<solid_energy_l_max_its",
      "The maximum allowed iterations in the linear solver of the solid energy equation.");

  params.addParamNamesToGroup("solid_energy_petsc_options solid_energy_petsc_options_iname "
                              "solid_energy_petsc_options_value solid_energy_absolute_tolerance "
                              "solid_energy_l_tol solid_energy_l_abs_tol solid_energy_l_max_its",
                              "Solid Energy Equation");

  /*
   * Parameters to control the solution of the phase equation
   */

  params.addRangeCheckedParam<Real>(
      "phase_equation_relaxation",
      0.9,
      "0.0<phase_equation_relaxation<=1.0",
      "The relaxation which should be used for the phase equation. (=1 for no relaxation, "
      "diagonal dominance will still be enforced)");

  params.addRangeCheckedParam<unsigned int>("MULES_iterations",
                                            1,
                                            "1<=MULES_iterations",
                                            "Number of MULES iterations to perform.");

  params.addRangeCheckedParam<unsigned int>(
      "n_alpha_subcycles",
      1,
      "n_alpha_subcycles >= 1",
      "Number of sub-cycles for the phase transport equation per PIMPLE iteration. "
      "The alpha equation is solved n_alpha_subcycles times with dt/n_alpha_subcycles.");

  params.addParam<MultiMooseEnum>("phase_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the phase equation");
  params.addParam<MultiMooseEnum>("phase_petsc_options_iname",
                                  Moose::PetscSupport::getCommonPetscKeys(),
                                  "Names of PETSc name/value pairs for the phase equation");
  params.addParam<std::vector<std::string>>(
      "phase_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "phase equation");

  params.addRangeCheckedParam<Real>(
      "phase_absolute_tolerance",
      1e-5,
      "0.0<phase_absolute_tolerance",
      "The absolute tolerance on the normalized residual of the energy equation.");

  params.addRangeCheckedParam<Real>("phase_l_tol",
                                    1e-5,
                                    "0.0<=phase_l_tol & phase_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the phase equation.");
  params.addRangeCheckedParam<Real>("phase_l_abs_tol",
                                    1e-10,
                                    "0.0<phase_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the phase equation.");
  params.addRangeCheckedParam<unsigned int>(
      "phase_l_max_its",
      10000,
      "0<phase_l_max_its",
      "The maximum allowed iterations in the linear solver of the phase equation.");

  params.addParamNamesToGroup(
      "phase_equation_relaxation MULES_iterations n_alpha_subcycles phase_petsc_options "
      "phase_petsc_options_iname "
      "phase_petsc_options_value phase_petsc_options_value phase_absolute_tolerance "
      "phase_l_tol phase_l_abs_tol phase_l_max_its",
      "Phase Equation");

  /*
   * Parameters to control the solution of each scalar advection system
   */
  params.addParam<std::vector<Real>>("passive_scalar_equation_relaxation",
                                     std::vector<Real>(),
                                     "The relaxation which should be used for the passive scalar "
                                     "equations. (=1 for no relaxation, "
                                     "diagonal dominance will still be enforced)");

  params.addParam<MultiMooseEnum>("passive_scalar_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the passive scalar equation(s)");
  params.addParam<MultiMooseEnum>(
      "passive_scalar_petsc_options_iname",
      Moose::PetscSupport::getCommonPetscKeys(),
      "Names of PETSc name/value pairs for the passive scalar equation(s)");
  params.addParam<std::vector<std::string>>(
      "passive_scalar_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "passive scalar equation(s)");
  params.addParam<std::vector<Real>>(
      "passive_scalar_absolute_tolerance",
      std::vector<Real>(),
      "The absolute tolerance(s) on the normalized residual(s) of the passive scalar equation(s).");
  params.addRangeCheckedParam<Real>("passive_scalar_l_tol",
                                    1e-5,
                                    "0.0<=passive_scalar_l_tol & passive_scalar_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the passive scalar equation(s).");
  params.addRangeCheckedParam<Real>("passive_scalar_l_abs_tol",
                                    1e-10,
                                    "0.0<passive_scalar_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the passive scalar equation(s).");
  params.addParam<unsigned int>(
      "passive_scalar_l_max_its",
      10000,
      "The maximum allowed iterations in the linear solver of the turbulence equation.");

  params.addParamNamesToGroup(
      "passive_scalar_systems passive_scalar_equation_relaxation passive_scalar_petsc_options "
      "passive_scalar_petsc_options_iname "
      "passive_scalar_petsc_options_value passive_scalar_petsc_options_value "
      "passive_scalar_absolute_tolerance "
      "passive_scalar_l_tol passive_scalar_l_abs_tol passive_scalar_l_max_its",
      "passive_scalar Equation");

  /*
   * Parameters to control the solution of each turbulence system
   */
  params.addParam<std::vector<Real>>("turbulence_equation_relaxation",
                                     std::vector<Real>(),
                                     "The relaxation which should be used for the turbulence "
                                     "equations. (=1 for no relaxation, "
                                     "diagonal dominance will still be enforced)");

  params.addParam<std::vector<Real>>("turbulence_field_relaxation",
                                     std::vector<Real>(),
                                     "The relaxation which should be used for the turbulence "
                                     "fields.");

  params.addParam<std::vector<Real>>(
      "turbulence_field_min_limit",
      std::vector<Real>(),
      "The lower limit imposed on turbulent quantities. The recommended value for robustness "
      "is 1e-8.");

  params.addParam<MultiMooseEnum>("turbulence_petsc_options",
                                  Moose::PetscSupport::getCommonPetscFlags(),
                                  "Singleton PETSc options for the turbulence equation(s)");
  params.addParam<MultiMooseEnum>("turbulence_petsc_options_iname",
                                  Moose::PetscSupport::getCommonPetscKeys(),
                                  "Names of PETSc name/value pairs for the turbulence equation(s)");
  params.addParam<std::vector<std::string>>(
      "turbulence_petsc_options_value",
      "Values of PETSc name/value pairs (must correspond with \"petsc_options_iname\" for the "
      "turbulence equation(s)");
  params.addParam<std::vector<Real>>(
      "turbulence_absolute_tolerance",
      std::vector<Real>(),
      "The absolute tolerance(s) on the normalized residual(s) of the turbulence equation(s).");
  params.addRangeCheckedParam<Real>("turbulence_l_tol",
                                    1e-5,
                                    "0.0<=turbulence_l_tol & turbulence_l_tol<1.0",
                                    "The relative tolerance on the normalized residual in the "
                                    "linear solver of the turbulence equation(s).");
  params.addRangeCheckedParam<Real>("turbulence_l_abs_tol",
                                    1e-10,
                                    "0.0<turbulence_l_abs_tol",
                                    "The absolute tolerance on the normalized residual in the "
                                    "linear solver of the turbulence equation(s).");
  params.addParam<unsigned int>(
      "turbulence_l_max_its",
      10000,
      "The maximum allowed iterations in the linear solver of the turbulence equation.");

  params.addParamNamesToGroup("turbulence_systems "
                              "turbulence_equation_relaxation "
                              "turbulence_field_relaxation "
                              "turbulence_field_min_limit "
                              "turbulence_petsc_options "
                              "turbulence_petsc_options_iname "
                              "turbulence_petsc_options_value turbulence_petsc_options_value "
                              "turbulence_absolute_tolerance "
                              "turbulence_l_tol turbulence_l_abs_tol turbulence_l_max_its",
                              "Turbulence Equations");

  /*
   * SIMPLE iteration control
   */

  params.addRangeCheckedParam<unsigned int>(
      "num_iterations",
      1000,
      "0<num_iterations",
      "The number of momentum-pressure-(other fields) iterations needed.");

  params.addParam<bool>("continue_on_max_its",
                        true,
                        "If solve should continue if maximum number of iterations is hit.");


  /*
   * Parameters to control interface sharpening
   */
  params.addParam<bool>("enforce_phase_sum", true, "Ensure that all phases add to 1.0 for a cell.");
  params.addParam<bool>("activate_interface_shapening", false, "Activate interface sharpening method.");
  MooseEnum sharpening_type("heaviside smooth", "heaviside");
  params.addParam<MooseEnum>("shapening_type", sharpening_type, "The interface sharpening method used.");
  params.addParam<Real>("smoothing_constant", 100.0, "The interface smoothing constant value.");
  
  params.addParamNamesToGroup("enforce_phase_sum "
                              "activate_interface_shapening "
                              "shapening_type "
                              "smoothing_constant ",
                              "Interface Sharpenning");
  return params;
}

SolveBaseMultiPhase::SolveBaseMultiPhase(Executioner & ex)
  : SolveObject(ex),
    UserObjectInterface(this),
    // Phases
    _number_of_phases(getParam<unsigned int>("number_of_phases")),
    // Momentum
    _momentum_system_names(getParam<std::vector<std::vector<std::string>>>("momentum_systems")),
    _momentum_l_abs_tol(getParam<Real>("momentum_l_abs_tol")),
    _momentum_equation_relaxation(getParam<Real>("momentum_equation_relaxation")),
    // Pressure
    _pressure_system_name(getParam<SolverSystemName>("pressure_system")),
    _pressure_l_abs_tol(getParam<Real>("pressure_l_abs_tol")),
    _pressure_variable_relaxation(getParam<Real>("pressure_variable_relaxation")),
    _pin_pressure(getParam<bool>("pin_pressure")),
    _pressure_pin_value(getParam<Real>("pressure_pin_value")),
    _pressure_pin_dof(libMesh::invalid_uint),
    // Fluid Energy
    _has_energy_system(isParamValid("energy_systems")),
    _energy_system_names(_has_energy_system ? &getParam<std::vector<SolverSystemName>>("energy_systems") : nullptr),
    _energy_equation_relaxation(getParam<Real>("energy_equation_relaxation")),
    _energy_l_abs_tol(getParam<Real>("energy_l_abs_tol")),
    // Solid Energy
    _has_solid_energy_system(_has_energy_system && isParamValid("solid_energy_system")),
    _solid_energy_system_names(isParamValid("solid_energy_system") ? getParam<SolverSystemName>("solid_energy_system") : SolverSystemName()),
    _solid_energy_equation_relaxation(getParam<Real>("solid_energy_equation_relaxation")),
    _solid_energy_l_abs_tol(getParam<Real>("solid_energy_l_abs_tol")),
    // Phase Transport - multiphase solver must have energy systems
    _phase_system_names(getParam<std::vector<SolverSystemName>>("phase_systems")),
    _number_of_solved_phases(_phase_system_names.size()),
    _phase_equation_relaxation(getParam<Real>("phase_equation_relaxation")),
    _MULES_iterations(getParam<unsigned int>("MULES_iterations")),
    _n_alpha_subcycles(getParam<unsigned int>("n_alpha_subcycles")),
    _phase_l_abs_tol(getParam<Real>("phase_l_abs_tol")),
    // Passive Scalars
    _passive_scalar_system_names(getParam<std::vector<SolverSystemName>>("passive_scalar_systems")),
    _has_passive_scalar_systems(!_passive_scalar_system_names.empty()),
    _passive_scalar_equation_relaxation(
        getParam<std::vector<Real>>("passive_scalar_equation_relaxation")),
    _passive_scalar_l_abs_tol(getParam<Real>("passive_scalar_l_abs_tol")),
    // Turbulence
    _turbulence_system_names(isParamValid("turbulence_systems") ? &getParam<std::vector<std::vector<std::string>>>("turbulence_systems") : nullptr),
    _has_turbulence_systems(isParamValid("turbulence_systems")),
    _turbulence_equation_relaxation(getParam<std::vector<Real>>("turbulence_equation_relaxation")),
    _turbulence_field_relaxation(getParam<std::vector<Real>>("turbulence_field_relaxation")),
    _turbulence_field_min_limit(getParam<std::vector<Real>>("turbulence_field_min_limit")),
    _turbulence_l_abs_tol(getParam<Real>("turbulence_l_abs_tol")),
    // Solve tolerances
    _momentum_absolute_tolerance(getParam<Real>("momentum_absolute_tolerance")),
    _pressure_absolute_tolerance(getParam<Real>("pressure_absolute_tolerance")),
    _energy_absolute_tolerance(getParam<Real>("energy_absolute_tolerance")),
    _solid_energy_absolute_tolerance(getParam<Real>("solid_energy_absolute_tolerance")),
    _phase_absolute_tolerance(getParam<Real>("phase_absolute_tolerance")),
    _passive_scalar_absolute_tolerance(
        getParam<std::vector<Real>>("passive_scalar_absolute_tolerance")),
    _turbulence_absolute_tolerance(getParam<std::vector<Real>>("turbulence_absolute_tolerance")),
    // Counting and output
    _num_iterations(getParam<unsigned int>("num_iterations")),
    _continue_on_max_its(getParam<bool>("continue_on_max_its")),
    _print_fields(getParam<bool>("print_fields")),
    // System numbering and storage
    _pressure_sys_number(_problem.linearSysNum(getParam<SolverSystemName>("pressure_system"))),
    _pressure_system(_problem.getLinearSystem(_pressure_sys_number)),
    _solid_energy_sys_number(
        _has_solid_energy_system
            ? _problem.linearSysNum(getParam<SolverSystemName>("solid_energy_system"))
            : libMesh::invalid_uint),
    _solid_energy_system(
        _has_solid_energy_system ? &_problem.getLinearSystem(_solid_energy_sys_number) : nullptr),
    // Interface sharpening
    _enforce_phase_sum(getParam<bool>("enforce_phase_sum")),
    _activate_interface_shapening(getParam<bool>("activate_interface_shapening")),
    _shapening_type(getParam<MooseEnum>("shapening_type")),
    _smoothing_constant(getParam<Real>("smoothing_constant"))
{
  // Checks errors and assemblies for the momentum system
  for(unsigned int i = 0; i < _number_of_phases; ++i)
    if (_momentum_system_names[i].size() != _problem.mesh().dimension())
      paramError("momentum_systems",
                "The number of momentum components should be equal to the number of "
                "spatial dimensions on the mesh.");

  const auto & momentum_petsc_options = getParam<MultiMooseEnum>("momentum_petsc_options");
  const auto & momentum_petsc_pair_options = getParam<MooseEnumItem, std::string>(
      "momentum_petsc_options_iname", "momentum_petsc_options_value");
  Moose::PetscSupport::addPetscFlagsToPetscOptions(
      momentum_petsc_options, "-", *this, _momentum_petsc_options);
  Moose::PetscSupport::addPetscPairsToPetscOptions(momentum_petsc_pair_options,
                                                   _problem.mesh().dimension(),
                                                   "-",
                                                   *this,
                                                   _momentum_petsc_options);

  _momentum_linear_control.real_valued_data["rel_tol"] = getParam<Real>("momentum_l_tol");
  _momentum_linear_control.real_valued_data["abs_tol"] = getParam<Real>("momentum_l_abs_tol");
  _momentum_linear_control.int_valued_data["max_its"] =
      getParam<unsigned int>("momentum_l_max_its");

  // Checks errors and assembly of the pressure system
  const auto & pressure_petsc_options = getParam<MultiMooseEnum>("pressure_petsc_options");
  const auto & pressure_petsc_pair_options = getParam<MooseEnumItem, std::string>(
      "pressure_petsc_options_iname", "pressure_petsc_options_value");
  Moose::PetscSupport::addPetscFlagsToPetscOptions(
      pressure_petsc_options, "-", *this, _pressure_petsc_options);
  Moose::PetscSupport::addPetscPairsToPetscOptions(pressure_petsc_pair_options,
                                                   _problem.mesh().dimension(),
                                                   "-",
                                                   *this,
                                                   _pressure_petsc_options);

  _pressure_linear_control.real_valued_data["rel_tol"] = getParam<Real>("pressure_l_tol");
  _pressure_linear_control.real_valued_data["abs_tol"] = getParam<Real>("pressure_l_abs_tol");
  _pressure_linear_control.int_valued_data["max_its"] =
      getParam<unsigned int>("pressure_l_max_its");

  // Checks errors and assembly of the energy system
  if (_has_energy_system)
  {
    const auto & energy_petsc_options = getParam<MultiMooseEnum>("energy_petsc_options");
    const auto & energy_petsc_pair_options = getParam<MooseEnumItem, std::string>(
        "energy_petsc_options_iname", "energy_petsc_options_value");
    Moose::PetscSupport::addPetscFlagsToPetscOptions(
        energy_petsc_options, "-", *this, _energy_petsc_options);
    Moose::PetscSupport::addPetscPairsToPetscOptions(
        energy_petsc_pair_options, _problem.mesh().dimension(), "-", *this, _energy_petsc_options);

    _energy_linear_control.real_valued_data["rel_tol"] = getParam<Real>("energy_l_tol");
    _energy_linear_control.real_valued_data["abs_tol"] = getParam<Real>("energy_l_abs_tol");
    _energy_linear_control.int_valued_data["max_its"] = getParam<unsigned int>("energy_l_max_its");
  }
  else
    checkDependentParameterError("energy_systems",
                                 {"energy_petsc_options",
                                  "energy_petsc_options_iname",
                                  "energy_petsc_options_value",
                                  "energy_l_tol",
                                  "energy_l_abs_tol",
                                  "energy_l_max_its",
                                  "energy_absolute_tolerance",
                                  "energy_equation_relaxation"},
                                 false);

  // Checks errors and assemblies of the solid energy system
  if (_has_solid_energy_system)
  {
    const auto & solid_energy_petsc_options =
        getParam<MultiMooseEnum>("solid_energy_petsc_options");
    const auto & solid_energy_petsc_pair_options = getParam<MooseEnumItem, std::string>(
        "solid_energy_petsc_options_iname", "solid_energy_petsc_options_value");
    Moose::PetscSupport::addPetscFlagsToPetscOptions(
        solid_energy_petsc_options, "-", *this, _solid_energy_petsc_options);
    Moose::PetscSupport::addPetscPairsToPetscOptions(solid_energy_petsc_pair_options,
                                                     _problem.mesh().dimension(),
                                                     "-",
                                                     *this,
                                                     _solid_energy_petsc_options);

    _solid_energy_linear_control.real_valued_data["rel_tol"] = getParam<Real>("solid_energy_l_tol");
    _solid_energy_linear_control.real_valued_data["abs_tol"] =
        getParam<Real>("solid_energy_l_abs_tol");
    _solid_energy_linear_control.int_valued_data["max_its"] =
        getParam<unsigned int>("solid_energy_l_max_its");
  }
  else
    checkDependentParameterError("solid_energy_system",
                                 {"solid_energy_petsc_options",
                                  "solid_energy_petsc_options_iname",
                                  "solid_energy_petsc_options_value",
                                  "solid_energy_l_tol",
                                  "solid_energy_l_abs_tol",
                                  "solid_energy_l_max_its",
                                  "solid_energy_absolute_tolerance",
                                  "solid_energy_equation_relaxation"},
                                 false);

  // Checks errors and assembly the phase system
  const auto & phase_petsc_options = getParam<MultiMooseEnum>("phase_petsc_options");
  const auto & phase_petsc_pair_options = getParam<MooseEnumItem, std::string>(
      "phase_petsc_options_iname", "phase_petsc_options_value");
  Moose::PetscSupport::addPetscFlagsToPetscOptions(
      phase_petsc_options, "-", *this, _phase_petsc_options);
  Moose::PetscSupport::addPetscPairsToPetscOptions(
      phase_petsc_pair_options, _problem.mesh().dimension(), "-", *this, _phase_petsc_options);

  _phase_linear_control.real_valued_data["rel_tol"] = getParam<Real>("phase_l_tol");
  _phase_linear_control.real_valued_data["abs_tol"] = getParam<Real>("phase_l_abs_tol");
  _phase_linear_control.int_valued_data["max_its"] = getParam<unsigned int>("phase_l_max_its");

  // Checks errors and assembly the passive scalar system
  if (_has_passive_scalar_systems)
  {
    if (_passive_scalar_system_names.size() != _passive_scalar_equation_relaxation.size())
      paramError("passive_scalar_equation_relaxation",
                 "The number of equation relaxation parameters does not match the number of "
                 "passive scalar equations!");
    if (_passive_scalar_system_names.size() != _passive_scalar_absolute_tolerance.size())
      paramError("passive_scalar_absolute_tolerance",
                 "The number of absolute tolerances does not match the number of "
                 "passive scalar equations!");

    const auto & passive_scalar_petsc_options =
        getParam<MultiMooseEnum>("passive_scalar_petsc_options");
    const auto & passive_scalar_petsc_pair_options = getParam<MooseEnumItem, std::string>(
        "passive_scalar_petsc_options_iname", "passive_scalar_petsc_options_value");
    Moose::PetscSupport::addPetscFlagsToPetscOptions(
        passive_scalar_petsc_options, "-", *this, _passive_scalar_petsc_options);
    Moose::PetscSupport::addPetscPairsToPetscOptions(passive_scalar_petsc_pair_options,
                                                     _problem.mesh().dimension(),
                                                     "-",
                                                     *this,
                                                     _passive_scalar_petsc_options);

    _passive_scalar_linear_control.real_valued_data["rel_tol"] =
        getParam<Real>("passive_scalar_l_tol");
    _passive_scalar_linear_control.real_valued_data["abs_tol"] =
        getParam<Real>("passive_scalar_l_abs_tol");
    _passive_scalar_linear_control.int_valued_data["max_its"] =
        getParam<unsigned int>("passive_scalar_l_max_its");
  }
  else
    checkDependentParameterError("passive_scalar_systems",
                                 {"passive_scalar_petsc_options",
                                  "passive_scalar_petsc_options_iname",
                                  "passive_scalar_petsc_options_value",
                                  "passive_scalar_l_tol",
                                  "passive_scalar_l_abs_tol",
                                  "passive_scalar_l_max_its",
                                  "passive_scalar_equation_relaxation",
                                  "passive_scalar_absolute_tolerance"},
                                 false);

  // Checks errors and assemblies the turbulence system
  if (_has_turbulence_systems)
  {
    for(unsigned int i = 0; i < _number_of_phases; ++i)
    {
      if ((*_turbulence_system_names)[i].size() != _turbulence_equation_relaxation.size())
        paramError("turbulence_equation_relaxation",
                  "The number of equation relaxation parameters does not match the number of "
                  "turbulence equations!");
      if ((*_turbulence_system_names)[i].size() != _turbulence_absolute_tolerance.size())
        paramError("turbulence_absolute_tolerance",
                  "The number of absolute tolerances does not match the number of "
                  "turbulence equations!");
    }

    if (_turbulence_field_min_limit.empty())
      // If no minimum bounds are given, initialize to default value 1e-8
      _turbulence_field_min_limit.resize((*_turbulence_system_names)[0].size(), 1e-8);

    // Assign turbulence field relaxation as 1.0 if not defined
    if (_turbulence_field_relaxation.empty())
      _turbulence_field_relaxation.resize((*_turbulence_system_names)[0].size(), 1.0);

    const auto & turbulence_petsc_options = getParam<MultiMooseEnum>("turbulence_petsc_options");
    const auto & turbulence_petsc_pair_options = getParam<MooseEnumItem, std::string>(
        "turbulence_petsc_options_iname", "turbulence_petsc_options_value");
    Moose::PetscSupport::addPetscFlagsToPetscOptions(
        turbulence_petsc_options, "-", *this, _turbulence_petsc_options);
    Moose::PetscSupport::addPetscPairsToPetscOptions(turbulence_petsc_pair_options,
                                                     _problem.mesh().dimension(),
                                                     "-",
                                                     *this,
                                                     _turbulence_petsc_options);

    _turbulence_linear_control.real_valued_data["rel_tol"] = getParam<Real>("turbulence_l_tol");
    _turbulence_linear_control.real_valued_data["abs_tol"] = getParam<Real>("turbulence_l_abs_tol");
    _turbulence_linear_control.int_valued_data["max_its"] =
        getParam<unsigned int>("turbulence_l_max_its");
  }
  else
    checkDependentParameterError("turbulence_systems",
                                 {"turbulence_petsc_options",
                                  "turbulence_petsc_options_iname",
                                  "turbulence_petsc_options_value",
                                  "turbulence_l_tol",
                                  "turbulence_l_abs_tol",
                                  "turbulence_l_max_its",
                                  "turbulence_equation_relaxation",
                                  "turbulence_field_relaxation",
                                  "turbulence_field_min_limit",
                                  "turbulence_absolute_tolerance"},
                                 false);

  // We fetch the systems and their numbers for the momentum equations.
  _momentum_system_numbers.resize(_number_of_phases);
  _momentum_systems.resize(_number_of_phases);
  for(unsigned int i = 0; i < _number_of_phases; ++i)
    for (auto system_i : index_range(_momentum_system_names[i]))
    {
      _momentum_system_numbers[i].push_back(_problem.linearSysNum(_momentum_system_names[i][system_i]));
      _momentum_systems[i].push_back(&_problem.getLinearSystem(_momentum_system_numbers[i][system_i]));
      _systems_to_solve.push_back(_momentum_systems[i].back());
    }

  _systems_to_solve.push_back(&_pressure_system);

  if (_has_energy_system)
    for(unsigned int i = 0; i < _number_of_phases; ++i)
    {
      _energy_system_numbers.push_back(_problem.linearSysNum((*_energy_system_names)[i]));
      _energy_systems.push_back(&_problem.getLinearSystem(_energy_system_numbers[i]));
      _systems_to_solve.push_back(_energy_systems[i]);
    }

  if (_has_solid_energy_system)
    _systems_to_solve.push_back(_solid_energy_system);

  // add systems for the phases transpored
  for(unsigned int i = 0; i < _number_of_solved_phases; ++i)
  {
    _phase_system_numbers.push_back(_problem.linearSysNum(_phase_system_names[i]));
    _phase_systems.push_back(&_problem.getLinearSystem(_phase_system_numbers[i]));
    _systems_to_solve.push_back(_phase_systems[i]);
  }

  // and for the turbulence surrogate equations
  if (_has_turbulence_systems)
  {
    _turbulence_system_numbers.resize(_number_of_phases);
    _turbulence_systems.resize(_number_of_phases);
    for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
    {
      for (auto system_i : index_range((*_turbulence_system_names)[phase_number]))
      {
        const auto linear_system_number = _problem.linearSysNum((*_turbulence_system_names)[phase_number][system_i]);
        _turbulence_system_numbers[phase_number].push_back(linear_system_number);
        _turbulence_systems[phase_number].push_back(&_problem.getLinearSystem(linear_system_number));
      }
    }
  }

  // and for the passive scalar equations
  if (_has_passive_scalar_systems)
    for (auto system_i : index_range(_passive_scalar_system_names))
    {
      _passive_scalar_system_numbers.push_back(
          _problem.linearSysNum(_passive_scalar_system_names[system_i]));
      _passive_scalar_systems.push_back(
          &_problem.getLinearSystem(_passive_scalar_system_numbers[system_i]));
      _systems_to_solve.push_back(_passive_scalar_systems.back());
    }

  // We disable the prefix here for the time being, the segregated solvers use a different approach
  // for setting the petsc parameters
  for (auto & system : _systems_to_solve)
    system->system().prefix_with_name(false);
}

void
SolveBaseMultiPhase::setupPressurePin()
{
  if (_pin_pressure)
    _pressure_pin_dof = NS::FV::findPointDoFID(_problem.getVariable(0, "pressure"),
                                               _problem.mesh(),
                                               getParam<Point>("pressure_pin_point"));
}

void
SolveBaseMultiPhase::checkDependentParameterError(const std::string & main_parameter,
                                              const std::vector<std::string> & dependent_parameters,
                                              const bool should_be_defined)
{
  for (const auto & param : dependent_parameters)
    if (parameters().isParamSetByUser(param) == !should_be_defined)
      paramError(param,
                 "This parameter should " + std::string(should_be_defined ? "" : "not") +
                     " be given by the user with the corresponding " + main_parameter +
                     " setting!");
}

void
SolveBaseMultiPhase::linkRhieChowUserObjects()
{
  for(unsigned int i = 0; i < _number_of_phases; ++i)
  {
    _rc_uo.push_back(const_cast<RhieChowMassFluxMultiPhase *>(&getUserObjectByName<RhieChowMassFluxMultiPhase>(getParam<std::vector<UserObjectName>>("rhie_chow_user_objects")[i])));
    _rc_uo[i]->linkMomentumPressureSystems(
        _momentum_systems[i], _pressure_system, _momentum_system_numbers[i]);

    // Initialize the face velocities in the RC object
    if (!_app.isRecovering())
    {
      _rc_uo[i]->initFaceMassFlux();
      // Initialize consistent mass flux from the initial alpha field so it's
      // available before the first momentum solve (required for alpha-first ordering)
      _rc_uo[i]->computeConsistentMassFlux();
    }
    _rc_uo[i]->initCouplingField();
  }
}

std::vector<std::pair<unsigned int, Real>>
SolveBaseMultiPhase::solveMomentumPredictor(const unsigned int phase_number)
{
  // Temporary storage for the (flux-normalized) residuals from
  // different momentum components
  std::vector<std::pair<unsigned int, Real>> its_normalized_residuals;

  LinearImplicitSystem & momentum_system_0 =
      libMesh::cast_ref<LinearImplicitSystem &>(_momentum_systems[phase_number][0]->system());

  PetscLinearSolver<Real> & momentum_solver =
      libMesh::cast_ref<PetscLinearSolver<Real> &>(*momentum_system_0.get_linear_solver());

  // Solve the momentum equations.
  // TO DO: These equations are VERY similar. If we can store the differences (things coming from
  // BCs for example) separately, it is enough to construct one matrix.
  for (const auto system_i : index_range(_momentum_systems[phase_number]))
  {
    _problem.setCurrentLinearSystem(_momentum_system_numbers[phase_number][system_i]);

    // We will need the right hand side and the solution of the next component
    LinearImplicitSystem & momentum_system =
        libMesh::cast_ref<LinearImplicitSystem &>(_momentum_systems[phase_number][system_i]->system());

    NumericVector<Number> & solution = *(momentum_system.solution);
    NumericVector<Number> & rhs = *(momentum_system.rhs);
    SparseMatrix<Number> & mmat = *(momentum_system.matrix);

    auto diff_diagonal = solution.zero_clone();

    // We assemble the matrix and the right hand side
    _problem.computeLinearSystemSys(momentum_system, mmat, rhs);

    // Still need to relax the right hand side with the same vector
    NS::FV::relaxMatrix(mmat, _momentum_equation_relaxation, *diff_diagonal);
    NS::FV::relaxRightHandSide(rhs, solution, *diff_diagonal);

    // The normalization factor depends on the right hand side so we need to recompute it for this
    // component
    Real norm_factor = NS::FV::computeNormalizationFactor(solution, mmat, rhs);

    // Very important, for deciding the convergence, we need the unpreconditioned
    // norms in the linear solve
    LibmeshPetscCall(KSPSetNormType(momentum_solver.ksp(), KSP_NORM_UNPRECONDITIONED));
    // Solve this component. We don't update the ghosted solution yet, that will come at the end
    // of the corrector step. Also setting the linear tolerances and maximum iteration counts.
    _momentum_linear_control.real_valued_data["abs_tol"] = _momentum_l_abs_tol * norm_factor;
    momentum_solver.set_solver_configuration(_momentum_linear_control);

    // We solve the equation
    auto its_resid_pair = momentum_solver.solve(mmat, mmat, solution, rhs);
    momentum_system.update();

    // We will reuse the preconditioner for every momentum system
    if (system_i == 0)
      momentum_solver.reuse_preconditioner(true);

    // Save the normalized residual
    its_normalized_residuals.push_back(
        std::make_pair(its_resid_pair.first, momentum_solver.get_initial_residual() / norm_factor));

    if (_print_fields)
    {
      _console << " solution after solve " << std::endl;
      solution.print();
      _console << " matrix when we solve " << std::endl;
      mmat.print();
      _console << " rhs when we solve " << std::endl;
      rhs.print();
      _console << " velocity solution component " << system_i << std::endl;
      solution.print();
      _console << "Norm factor " << norm_factor << std::endl;
      _console << Moose::stringify(momentum_solver.get_initial_residual()) << std::endl;
    }

    // Printing residuals
    _console << " Momentum equation - " << "Phase " << phase_number << ": "
             << (_momentum_systems.size() > 1
                     ? std::string(" Component ") + std::to_string(system_i + 1) + std::string(" ")
                     : std::string(" "))
             << COLOR_GREEN << its_normalized_residuals[system_i].second << COLOR_DEFAULT
             << " Linear its: " << its_normalized_residuals[system_i].first << std::endl;
  }

  for (const auto system_i : index_range(_momentum_systems[phase_number]))
  {
    LinearImplicitSystem & momentum_system =
        libMesh::cast_ref<LinearImplicitSystem &>(_momentum_systems[phase_number][system_i]->system());
    _momentum_systems[phase_number][system_i]->setSolution(*(momentum_system.current_local_solution));
    _momentum_systems[phase_number][system_i]->copyPreviousNonlinearSolutions();
  }

  // We reset this to ensure the preconditioner is recomputed new time we go to the momentum
  // predictor
  momentum_solver.reuse_preconditioner(false);

  return its_normalized_residuals;
}

std::pair<unsigned int, Real>
SolveBaseMultiPhase::solvePressureCorrector()
{
  _problem.setCurrentLinearSystem(_pressure_sys_number);

  // We will need some members from the linear system
  LinearImplicitSystem & pressure_system =
      libMesh::cast_ref<LinearImplicitSystem &>(_pressure_system.system());

  // We will need the solution, the right hand side and the matrix
  NumericVector<Number> & current_local_solution = *(pressure_system.current_local_solution);
  NumericVector<Number> & solution = *(pressure_system.solution);
  SparseMatrix<Number> & mmat = *(pressure_system.matrix);
  NumericVector<Number> & rhs = *(pressure_system.rhs);

  // Fetch the linear solver from the system
  PetscLinearSolver<Real> & pressure_solver =
      libMesh::cast_ref<PetscLinearSolver<Real> &>(*pressure_system.get_linear_solver());

  _problem.computeLinearSystemSys(pressure_system, mmat, rhs, false);

  if (_print_fields)
  {
    _console << "Pressure matrix" << std::endl;
    mmat.print();
  }

  // We compute the normalization factors based on the fluxes
  Real norm_factor = NS::FV::computeNormalizationFactor(solution, mmat, rhs);

  // We need the non-preconditioned norm to be consistent with the norm factor
  LibmeshPetscCall(KSPSetNormType(pressure_solver.ksp(), KSP_NORM_UNPRECONDITIONED));

  // Setting the linear tolerances and maximum iteration counts
  _pressure_linear_control.real_valued_data["abs_tol"] = _pressure_l_abs_tol * norm_factor;
  pressure_solver.set_solver_configuration(_pressure_linear_control);

  if (_pin_pressure)
    NS::FV::constrainSystem(mmat, rhs, _pressure_pin_value, _pressure_pin_dof);
  pressure_system.update();

  auto its_res_pair = pressure_solver.solve(mmat, mmat, solution, rhs);
  pressure_system.update();

  if (_print_fields)
  {
    _console << " rhs when we solve pressure " << std::endl;
    rhs.print();
    _console << " Pressure " << std::endl;
    solution.print();
    _console << "Norm factor " << norm_factor << std::endl;
  }

  _pressure_system.setSolution(current_local_solution);

  const auto residuals =
      std::make_pair(its_res_pair.first, pressure_solver.get_initial_residual() / norm_factor);

  _console << " Pressure equation: " << COLOR_GREEN << residuals.second << COLOR_DEFAULT
           << " Linear its: " << residuals.first << std::endl;

  return residuals;
}

std::pair<unsigned int, Real>
SolveBaseMultiPhase::solveSolidEnergy()
{
  _problem.setCurrentLinearSystem(_solid_energy_sys_number);

  // We will need some members from the linear system
  LinearImplicitSystem & system =
      libMesh::cast_ref<LinearImplicitSystem &>(_solid_energy_system->system());

  // We will need the solution, the right hand side and the matrix
  NumericVector<Number> & current_local_solution = *(system.current_local_solution);
  NumericVector<Number> & solution = *(system.solution);
  SparseMatrix<Number> & mmat = *(system.matrix);
  NumericVector<Number> & rhs = *(system.rhs);

  // Fetch the linear solver from the system
  PetscLinearSolver<Real> & solver =
      libMesh::cast_ref<PetscLinearSolver<Real> &>(*system.get_linear_solver());

  _problem.computeLinearSystemSys(system, mmat, rhs, false);

  if (_print_fields)
  {
    _console << "Solid energy matrix" << std::endl;
    mmat.print();
  }

  // We compute the normalization factors based on the fluxes
  Real norm_factor = NS::FV::computeNormalizationFactor(solution, mmat, rhs);

  // We need the non-preconditioned norm to be consistent with the norm factor
  LibmeshPetscCall(KSPSetNormType(solver.ksp(), KSP_NORM_UNPRECONDITIONED));

  // Setting the linear tolerances and maximum iteration counts
  _solid_energy_linear_control.real_valued_data["abs_tol"] = _solid_energy_l_abs_tol * norm_factor;
  solver.set_solver_configuration(_solid_energy_linear_control);

  auto its_res_pair = solver.solve(mmat, mmat, solution, rhs);
  system.update();

  if (_print_fields)
  {
    _console << " rhs when we solve solid energy " << std::endl;
    rhs.print();
    _console << " Solid energy " << std::endl;
    solution.print();
    _console << "Norm factor " << norm_factor << std::endl;
  }

  _solid_energy_system->setSolution(current_local_solution);

  const auto residuals =
      std::make_pair(its_res_pair.first, solver.get_initial_residual() / norm_factor);

  _console << " Solid energy equation: " << COLOR_GREEN << residuals.second << COLOR_DEFAULT
           << " Linear its: " << residuals.first << std::endl;

  return residuals;
}

std::pair<unsigned int, Real>
SolveBaseMultiPhase::correctVelocities(const bool subtract_updated_pressure,
                             const bool recompute_face_mass_flux,
                             const SolverParams & solver_params)
{
  // Compute the coupling fields between the momentum and pressure equations.
  // The first argument makes sure the pressure gradient is staged at the first
  // iteration
  for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
    _rc_uo[phase_number]->computeHbyA(subtract_updated_pressure, _print_fields);

  // We set the preconditioner/controllable parameters for the pressure equations through
  // petsc options. Linear tolerances will be overridden within the solver.
  Moose::PetscSupport::petscSetOptions(_pressure_petsc_options, solver_params);

  // Solve the pressure corrector
  const auto residuals = solvePressureCorrector();

  // Compute the face velocity which is used in the advection terms. In certain
  // segregated solver algorithms (like PISO) this is only done on the last iteration.
  if (recompute_face_mass_flux)
  {
    for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
      _rc_uo[phase_number]->computeFaceMassFlux();
  }

  auto & pressure_current_solution = *(_pressure_system.system().current_local_solution.get());
  auto & pressure_old_solution = *(_pressure_system.solutionPreviousNewton());

  // Relax the pressure update for the next momentum predictor
  NS::FV::relaxSolutionUpdate(
      pressure_current_solution, pressure_old_solution, _pressure_variable_relaxation);

  // Overwrite old solution
  pressure_old_solution = pressure_current_solution;
  _pressure_system.setSolution(pressure_current_solution);

  // We recompute the updated pressure gradient
  _pressure_system.computeGradients();

  // Reconstruct the cell velocity as well to accelerate convergence
  for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
    _rc_uo[phase_number]->computeCellVelocity();

  return residuals;
}

std::pair<unsigned int, Real>
SolveBaseMultiPhase::solveAdvectedSystem(const unsigned int system_num,
                                                   LinearSystem & system,
                                                   const Real relaxation_factor,
                                                   libMesh::SolverConfiguration & solver_config,
                                                   const Real absolute_tol,
                                                   const bool relax_fields,
                                                   const Real field_relaxation)
{
  _problem.setCurrentLinearSystem(system_num);

  // We will need some members from the implicit linear system
  LinearImplicitSystem & li_system = libMesh::cast_ref<LinearImplicitSystem &>(system.system());

  // We will need the solution, the right hand side and the matrix
  NumericVector<Number> & current_local_solution = *(li_system.current_local_solution);
  NumericVector<Number> & solution = *(li_system.solution);
  SparseMatrix<Number> & mmat = *(li_system.matrix);
  NumericVector<Number> & rhs = *(li_system.rhs);

  // We need a vector that stores the (diagonal_relaxed-original_diagonal) vector
  auto diff_diagonal = solution.zero_clone();

  // Fetch the linear solver from the system
  PetscLinearSolver<Real> & linear_solver =
      libMesh::cast_ref<PetscLinearSolver<Real> &>(*li_system.get_linear_solver());

  _problem.computeLinearSystemSys(li_system, mmat, rhs, true);

  // Go and relax the system matrix and the right hand side
  NS::FV::relaxMatrix(mmat, relaxation_factor, *diff_diagonal);
  NS::FV::relaxRightHandSide(rhs, solution, *diff_diagonal);

  if (_print_fields)
  {
    _console << system.name() << " system matrix" << std::endl;
    mmat.print();
  }

  // We compute the normalization factors based on the fluxes
  Real norm_factor = NS::FV::computeNormalizationFactor(solution, mmat, rhs);

  // We need the non-preconditioned norm to be consistent with the norm factor
  LibmeshPetscCall(KSPSetNormType(linear_solver.ksp(), KSP_NORM_UNPRECONDITIONED));

  // Setting the linear tolerances and maximum iteration counts
  solver_config.real_valued_data["abs_tol"] = absolute_tol * norm_factor;
  linear_solver.set_solver_configuration(solver_config);

  // Solve the system and update current local solution
  auto its_res_pair = linear_solver.solve(mmat, mmat, solution, rhs);
  li_system.update();

  if (_print_fields)
  {
    _console << " rhs when we solve " << system.name() << std::endl;
    rhs.print();
    _console << system.name() << " solution " << std::endl;
    solution.print();
    _console << " Norm factor " << norm_factor << std::endl;
  }

  // Relax the field update for the next momentum predictor
  if (relax_fields)
  {
    const auto & old_local_solution = *(system.solutionPreviousNewton());
    NS::FV::relaxSolutionUpdate(current_local_solution, old_local_solution, field_relaxation);
  }

  system.setSolution(current_local_solution);

  const auto residuals =
      std::make_pair(its_res_pair.first, linear_solver.get_initial_residual() / norm_factor);

  _console << " Advected system: " << system.name() << " " << COLOR_GREEN << residuals.second
           << COLOR_DEFAULT << " Linear its: " << residuals.first << std::endl;

  return residuals;
}

bool
SolveBaseMultiPhase::solve()
{
  // Do not solve if problem is set not to
  if (!_problem.shouldSolve())
    return true;

  // ------------------------------------------------------------------
  //  Look up alpha advection kernels (once) for semi-implicit MULES
  // ------------------------------------------------------------------
  if (_alpha_advection_kernels.empty() && _number_of_solved_phases > 0)
  {
    for (unsigned int phase_number = 0; phase_number < _number_of_solved_phases; ++phase_number)
    {
      LinearFVMultiPhaseFractionAdvection * found_kernel = nullptr;
      std::vector<LinearFVFluxKernel *> kernels;
      _problem.theWarehouse()
          .query()
          .template condition<AttribSysNum>(_phase_system_numbers[phase_number])
          .template condition<AttribSystem>("LinearFVFluxKernel")
          .template condition<AttribThread>(0)
          .queryInto(kernels);

      for (auto * kernel : kernels)
      {
        auto * alpha_kernel =
            dynamic_cast<LinearFVMultiPhaseFractionAdvection *>(kernel);
        if (alpha_kernel && alpha_kernel->isSemiImplicitMULES())
        {
          found_kernel = alpha_kernel;
          break;
        }
      }
      _alpha_advection_kernels.push_back(found_kernel);
    }
  }

  // ------------------------------------------------------------------
  //  Capture old face flux for ddtCorr (once per timestep, before PIMPLE)
  // ------------------------------------------------------------------
  for (unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
    _rc_uo[phase_number]->captureOldFlux();

  // ------------------------------------------------------------------
  //  Helper counts of equations per phase
  // ------------------------------------------------------------------
  const unsigned int n_vel  = _momentum_systems.front().size();
  const unsigned int n_turb = _has_turbulence_systems
                                ? (*_turbulence_system_names)[0].size()
                                : 0;

  // ------------------------------------------------------------------
  //  Residual and tolerance vectors – sized *exactly* once
  // ------------------------------------------------------------------
  unsigned int no_systems =
        n_vel * _number_of_phases                       // momentum
      + 1                                               // pressure
      + (_has_energy_system       ? _number_of_phases : 0)
      + (_has_solid_energy_system ? 1                 : 0)
      + _number_of_solved_phases                        // phase fractions
      + n_turb * _number_of_phases;                    // turbulence

  std::vector<std::pair<unsigned int, Real>> ns_residuals(no_systems, std::make_pair(0u, 1.0));
  std::vector<Real> ns_abs_tols;
  ns_abs_tols.reserve(no_systems);

  // phases (solved first in alpha-first order)
  for (unsigned int i = 0; i < _number_of_solved_phases; ++i)
    ns_abs_tols.push_back(_phase_absolute_tolerance);

  // momentum tolerances
  for (unsigned int p = 0; p < _number_of_phases; ++p)
    for (unsigned int c = 0; c < n_vel; ++c)
      ns_abs_tols.push_back(_momentum_absolute_tolerance);

  // pressure
  ns_abs_tols.push_back(_pressure_absolute_tolerance);

  // energy (fluid)
  if (_has_energy_system)
    for (unsigned int p = 0; p < _number_of_phases; ++p)
      ns_abs_tols.push_back(_energy_absolute_tolerance);

  // energy (solid)
  if (_has_solid_energy_system)
    ns_abs_tols.push_back(_solid_energy_absolute_tolerance);

  // turbulence
  if (_has_turbulence_systems)
    for (unsigned int p = 0; p < _number_of_phases; ++p)
      for (unsigned int eq = 0; eq < n_turb; ++eq)
        ns_abs_tols.push_back(_turbulence_absolute_tolerance[eq]);

  // ------------------------------------------------------------------
  //  Constant solver parameters
  // ------------------------------------------------------------------
  SolverParams solver_params;
  solver_params._type = Moose::SolveType::ST_LINEAR;
  solver_params._line_search = Moose::LineSearchType::LS_NONE;

  // ------------------------------------------------------------------
  //  SIMPLE / PIMPLE outer loop
  // ------------------------------------------------------------------
  unsigned int simple_iteration_counter = 0;
  bool converged = false;

  // Loop until converged or hit the maximum allowed iteration number
  while (simple_iteration_counter < _num_iterations && !converged)
  {
    simple_iteration_counter++;

    // We set the preconditioner/controllable parameters through petsc options. Linear
    // tolerances will be overridden within the solver. In case of a segregated momentum
    // solver, we assume that every velocity component uses the same preconditioner
    Moose::PetscSupport::petscSetOptions(_momentum_petsc_options, solver_params);

    // Initialize pressure gradients, after this we just reuse the last ones from each
    // iteration
    if (simple_iteration_counter == 1)
      _pressure_system.computeGradients();

    _console << "Iteration " << simple_iteration_counter << " Initial residual norms:" << std::endl;

    unsigned int residual_counter = 0;

    // ---------------------------------------------------------------
    // 1. Phase transport with optional MULES sub-iterations
    //    and alpha sub-cycling. Alpha-first: solve before momentum
    //    so mixture properties are current for momentum/pressure.
    // ---------------------------------------------------------------
    const auto residual_counter_base = residual_counter;

    // Save original dt and apply sub-cycling reduction
    const Real dt_global = _problem.dt();
    if (_n_alpha_subcycles > 1)
      _problem.dt() = dt_global / static_cast<Real>(_n_alpha_subcycles);

    for (unsigned int subcycle = 0; subcycle < _n_alpha_subcycles; ++subcycle)
    {
      if (_n_alpha_subcycles > 1)
        _console << COLOR_CYAN << " -- ALPHA SUB-CYCLE: " << subcycle + 1
                 << " / " << _n_alpha_subcycles << COLOR_DEFAULT << std::endl;

      for(unsigned int MULES_iteration = 0; MULES_iteration < _MULES_iterations; ++MULES_iteration)
      {
        if(_MULES_iterations > 1)
          _console << COLOR_CYAN << " -- MULES ITERATION: " << MULES_iteration << std::endl;

        for(unsigned int phase_number = 0; phase_number < _number_of_solved_phases; ++phase_number)
        {
          // We set the preconditioner/controllable parameters through petsc options. Linear
          // tolerances will be overridden within the solver.
          Moose::PetscSupport::petscSetOptions(_phase_petsc_options, solver_params);
          ns_residuals[residual_counter_base + phase_number] =
              solveAdvectedSystem(_phase_system_numbers[phase_number],
                                  *_phase_systems[phase_number],
                                  _phase_equation_relaxation,
                                  _phase_linear_control,
                                  _phase_l_abs_tol);

          // Apply semi-implicit MULES correction (second pass: bounded anti-diffusive correction)
          // Applied on the last MULES iteration of every sub-cycle (matching OpenFOAM interFoam).
          if (MULES_iteration == _MULES_iterations - 1 &&
              phase_number < _alpha_advection_kernels.size() &&
              _alpha_advection_kernels[phase_number])
          {
            _alpha_advection_kernels[phase_number]->applyMULESCorrection(
                *_phase_systems[phase_number]);
            // Clamp alpha to [0,1] immediately after correction — the per-face MULES
            // limiter doesn't guarantee per-cell boundedness when corrections accumulate.
            LinearImplicitSystem & li_sys =
                libMesh::cast_ref<LinearImplicitSystem &>(_phase_systems[phase_number]->system());
            NS::FV::limitSolutionUpdate(*(li_sys.solution), 0.0, 1.0);
          }

          // Update residual counter only on the last MULES iteration of the last sub-cycle
          if (MULES_iteration == _MULES_iterations - 1 && subcycle == _n_alpha_subcycles - 1)
            residual_counter++;
        }

        // Bound individual phases
        for(unsigned int phase_number = 0; phase_number < _number_of_solved_phases; ++phase_number)
        {
          LinearImplicitSystem & li_system =
              libMesh::cast_ref<LinearImplicitSystem &>(_phase_systems[phase_number]->system());
          NumericVector<Number> & current_solution = *(li_system.solution);
          NS::FV::limitSolutionUpdate(current_solution, 0.0, 1.0);
        }
      }

      // After each sub-cycle, copy current solution to "old" for the next sub-cycle
      // so the time derivative kernel uses the just-computed alpha as alpha_old
      if (subcycle < _n_alpha_subcycles - 1)
      {
        for(unsigned int phase_number = 0; phase_number < _number_of_solved_phases; ++phase_number)
        {
          LinearImplicitSystem & li_system =
              libMesh::cast_ref<LinearImplicitSystem &>(_phase_systems[phase_number]->system());
          _phase_systems[phase_number]->setSolution(*(li_system.solution));
          _phase_systems[phase_number]->copyPreviousNonlinearSolutions();
        }
      }
    }

    // Restore global dt
    if (_n_alpha_subcycles > 1)
      _problem.dt() = dt_global;

    // Apply interface sharpening (inside PIMPLE loop after phase solve)
    if(_activate_interface_shapening)
    {
      for(unsigned int phase_number = 0; phase_number < _number_of_solved_phases; ++phase_number)
      {
        LinearImplicitSystem & li_system =
            libMesh::cast_ref<LinearImplicitSystem &>(_phase_systems[phase_number]->system());
        NumericVector<Number> & current_solution = *(li_system.solution);
        NS::FV::sharpenPhaseField(current_solution,
                                  *(_rc_uo[0]->getCellVolumes()),
                                  _shapening_type,
                                  _smoothing_constant);
      }
    }

    // Enforce phase sum constraint
    if (_number_of_phases == _number_of_solved_phases && _enforce_phase_sum)
    {
      std::vector<NumericVector<Number> *> phase_solutions;
      phase_solutions.reserve(_number_of_solved_phases);
      for(unsigned int phase_number = 0; phase_number < _number_of_solved_phases; ++phase_number)
      {
        LinearImplicitSystem & li_system =
            libMesh::cast_ref<LinearImplicitSystem &>(_phase_systems[phase_number]->system());
        NumericVector<Number> & current_solution = *(li_system.solution);
        phase_solutions.push_back(&current_solution);
      }
      NS::FV::constrainPhaseUpdate(phase_solutions);
    }

    // Compute consistent mass flux for momentum equation (rhoPhi from MULES-limited face alpha)
    for (unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
      _rc_uo[phase_number]->computeConsistentMassFlux();

    // Update material properties (rho_mix, mu_mix) from new alpha
    _problem.execute(EXEC_NONLINEAR);

    // ---------------------------------------------------------------
    // 2. Momentum predictor
    // ---------------------------------------------------------------
    // Solve the momentum predictor step
    for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
    {
      auto momentum_residual = solveMomentumPredictor(phase_number);
      for (const auto system_i : index_range(momentum_residual))
      {
        ns_residuals[residual_counter] = momentum_residual[system_i];
        residual_counter++;
      }
    }

    // ---------------------------------------------------------------
    // 3. Pressure corrector (and cell/face velocity update)
    // ---------------------------------------------------------------
    // Now we correct the velocity, this function depends on the method, it differs for
    // SIMPLE/PIMPLE, this returns the pressure errors
    ns_residuals[residual_counter] = correctVelocities(true, true, solver_params);
    residual_counter++;

    // ---------------------------------------------------------------
    // 4. Fluid energy (per phase)
    // ---------------------------------------------------------------
    // If we have an energy equation, solve it here. We assume the material properties in the
    // Navier-Stokes equations depend on temperature, therefore we can not solve for temperature
    // outside of the velocity-pressure loop
    if (_has_energy_system)
    {
      for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
      {
        // We set the preconditioner/controllable parameters through petsc options. Linear
        // tolerances will be overridden within the solver.
        Moose::PetscSupport::petscSetOptions(_energy_petsc_options, solver_params);
        ns_residuals[residual_counter] =
            solveAdvectedSystem(_energy_system_numbers[phase_number],
                                *_energy_systems[phase_number],
                                _energy_equation_relaxation,
                                _energy_linear_control,
                                _energy_l_abs_tol);
        residual_counter++;
      }
    }

    // ---------------------------------------------------------------
    // 5. Solid energy
    // ---------------------------------------------------------------
    if (_has_solid_energy_system)
    {
      // We set the preconditioner/controllable parameters through petsc options. Linear
      // tolerances will be overridden within the solver.
      Moose::PetscSupport::petscSetOptions(_solid_energy_petsc_options, solver_params);
      ns_residuals[residual_counter] = solveSolidEnergy();
      residual_counter++;
    }

    // ---------------------------------------------------------------
    // 6. Turbulence surrogate equations
    // ---------------------------------------------------------------
    // If we have turbulence equations, solve them here.
    // The turbulent viscosity depends on the value of the turbulence surrogate variables
    if (_has_turbulence_systems)
    {
      // We set the preconditioner/controllable parameters through petsc options. Linear
      // tolerances will be overridden within the solver.
      Moose::PetscSupport::petscSetOptions(_turbulence_petsc_options, solver_params);

      for(unsigned int phase_number = 0; phase_number < _number_of_phases; ++phase_number)
      {
        for (const auto i : index_range((*_turbulence_system_names)[phase_number]))
        {
            ns_residuals[residual_counter] =
                solveAdvectedSystem(_turbulence_system_numbers[phase_number][i],
                                    *_turbulence_systems[phase_number][i],
                                    _turbulence_equation_relaxation[i],
                                    _turbulence_linear_control,
                                    _turbulence_l_abs_tol,
                                    true,
                                    _turbulence_field_relaxation[i]);
            residual_counter++;

            // Limiting turbulence solution
            LinearImplicitSystem & li_system =
                libMesh::cast_ref<LinearImplicitSystem &>(_turbulence_systems[phase_number][i]->system());
            NumericVector<Number> & current_solution = *(li_system.solution);
            NS::FV::limitSolutionUpdate(current_solution, _turbulence_field_min_limit[i]);

            // Overwrite old solution
            _turbulence_systems[phase_number][i]->setSolution(current_solution);
          }
      }
    }

    // ---------------------------------------------------------------
    // 7. Check convergence of the flow block
    // ---------------------------------------------------------------
    converged = NS::FV::converged(ns_residuals, ns_abs_tols);
  }

  // ------------------------------------------------------------------
  // 10. Passive scalars (outside main loop)
  // ------------------------------------------------------------------
  // If we have passive scalar equations, solve them here. We assume the material properties in the
  // Navier-Stokes equations do not depend on passive scalars, as they are passive, therefore we
  // solve outside of the velocity-pressure loop
  if (_has_passive_scalar_systems && (converged || _continue_on_max_its))
  {
    // The reason why we need more than one iteration is due to the matrix relaxation
    // which can be used to stabilize the equations
    bool passive_scalar_converged = false;
    unsigned int ps_iteration_counter = 0;

    _console << "Passive scalar iteration " << ps_iteration_counter
             << " Initial residual norms:" << std::endl;

    while (ps_iteration_counter < _num_iterations && !passive_scalar_converged)
    {
      ps_iteration_counter++;
      std::vector<std::pair<unsigned int, Real>> scalar_residuals(
          _passive_scalar_system_names.size(), std::make_pair(0, 1.0));
      std::vector<Real> scalar_abs_tols;
      for (const auto scalar_tol : _passive_scalar_absolute_tolerance)
        scalar_abs_tols.push_back(scalar_tol);

      // We set the preconditioner/controllable parameters through petsc options. Linear
      // tolerances will be overridden within the solver.
      Moose::PetscSupport::petscSetOptions(_passive_scalar_petsc_options, solver_params);
      for (const auto i : index_range(_passive_scalar_system_names))
        scalar_residuals[i] = solveAdvectedSystem(_passive_scalar_system_numbers[i],
                                                  *_passive_scalar_systems[i],
                                                  _passive_scalar_equation_relaxation[i],
                                                  _passive_scalar_linear_control,
                                                  _passive_scalar_l_abs_tol);

      passive_scalar_converged = NS::FV::converged(scalar_residuals, scalar_abs_tols);
    }

    // Both flow and scalars must converge
    converged = passive_scalar_converged && converged;
  }

  converged = _continue_on_max_its ? true : converged;

  return converged;
}

