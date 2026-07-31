//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LowMachImplicitMassFlux.h"

#include "LinearSystem.h"
#include "LinearFVLowMachMassAdvection.h"
#include "MooseLinearVariableFV.h"
#include "TheWarehouse.h"

#include <algorithm>
#include <cmath>
#include <limits>

registerMooseObject("NavierStokesApp", LowMachImplicitMassFlux);

InputParameters
LowMachImplicitMassFlux::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params.addClassDescription(
      "Publishes the converged temporary-density mass flux from the same discrete face operator "
      "used by the implicit low-Mach mass equation.");
  params.addRequiredParam<SolverSystemName>("system", "Temporary-density linear system.");
  params.addRequiredParam<std::string>("mass_advection_kernel",
                                       "Name of the LinearFVLowMachMassAdvection kernel.");
  params.addRequiredParam<MooseFunctorName>(
      "temporary_density", "Temporary-density variable transported by the mass equation.");
  params.addRequiredParam<MooseFunctorName>(
      "eos_density", "Thermodynamic EOS density used to synchronize the accepted state.");
  params.addParam<MooseFunctorName>(
      "mass_flux_name", "low_mach_mass_flux", "Name of the published face mass-flux functor.");
  params.addRangeCheckedParam<Real>(
      "minimum_density", 0.0, "minimum_density>=0", "Minimum physically admissible density.");
  params.addParam<Real>("maximum_density",
                        std::numeric_limits<Real>::max(),
                        "Maximum physically admissible density.");
  params.addRangeCheckedParam<Real>(
      "density_tolerance",
      1e-12,
      "density_tolerance>=0",
      "Tolerance outside the admissible density bounds before reporting an error.");
  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_NONE);
  exec_enum = {EXEC_NONE};
  return params;
}

LowMachImplicitMassFlux::LowMachImplicitMassFlux(const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    _system_name(getParam<SolverSystemName>("system")),
    _mass_advection_kernel_name(getParam<std::string>("mass_advection_kernel")),
    _density_variable_name(getParam<MooseFunctorName>("temporary_density")),
    _minimum_density(getParam<Real>("minimum_density")),
    _maximum_density(getParam<Real>("maximum_density")),
    _density_tolerance(getParam<Real>("density_tolerance")),
    _eos_density(getFunctor<Real>("eos_density")),
    _mass_flux_name(getParam<MooseFunctorName>("mass_flux_name")),
    _mass_flux(_fe_problem.mesh(), _mass_flux_name)
{
  if (_maximum_density < _minimum_density)
    paramError("maximum_density", "maximum_density must be no smaller than minimum_density.");

  for (const auto tid : make_range(libMesh::n_threads()))
    UserObject::_subproblem.addFunctor(_mass_flux_name, _mass_flux, tid);
}

void
LowMachImplicitMassFlux::initialSetup()
{
  linkMassAdvectionKernel();
  initializeFluxStorage();
}

void
LowMachImplicitMassFlux::meshChanged()
{
  linkMassAdvectionKernel();
  initializeFluxStorage();
}

void
LowMachImplicitMassFlux::linkMassAdvectionKernel()
{
  const auto system_number = _fe_problem.linearSysNum(_system_name);
  std::vector<LinearFVFluxKernel *> kernels;
  _fe_problem.theWarehouse()
      .query()
      .condition<AttribThread>(0)
      .condition<AttribSysNum>(system_number)
      .condition<AttribSystem>("LinearFVFluxKernel")
      .condition<AttribName>(_mass_advection_kernel_name)
      .queryInto(kernels);

  if (kernels.size() != 1)
    paramError("mass_advection_kernel",
               "Expected exactly one linear FV flux kernel named '",
               _mass_advection_kernel_name,
               "' in system '",
               _system_name,
               "'.");

  _mass_advection_kernel = dynamic_cast<LinearFVLowMachMassAdvection *>(kernels.front());
  if (!_mass_advection_kernel)
    paramError("mass_advection_kernel",
               "Kernel '",
               _mass_advection_kernel_name,
               "' must be a LinearFVLowMachMassAdvection.");
  if (_mass_advection_kernel->variable().name() != _density_variable_name)
    paramError("temporary_density",
               "The named temporary density must be the variable of mass-advection kernel '",
               _mass_advection_kernel_name,
               "'.");
}

void
LowMachImplicitMassFlux::initializeFluxStorage()
{
  _mass_flux.clear();
  for (const auto * face_info : _fe_problem.mesh().faceInfo())
    _mass_flux[face_info->id()] = 0.0;
}

void
LowMachImplicitMassFlux::execute()
{
  publishConvergedMassFlux();
}

void
LowMachImplicitMassFlux::publishConvergedMassFlux()
{
  if (!_mass_advection_kernel)
    mooseError(name(), ": mass-advection kernel has not been linked.");

  const auto & variable = _mass_advection_kernel->variable();
  const auto state = Moose::currentState();
  Real local_minimum = std::numeric_limits<Real>::max();
  Real local_maximum = std::numeric_limits<Real>::lowest();
  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
    if (_mass_advection_kernel->hasBlocks(elem_info->subdomain_id()))
    {
      const Real density = variable.getElemValue(*elem_info, state);
      local_minimum = std::min(local_minimum, density);
      local_maximum = std::max(local_maximum, density);
    }

  Real global_minimum = local_minimum;
  Real global_maximum = local_maximum;
  _communicator.min(global_minimum);
  _communicator.max(global_maximum);
  if (global_minimum < _minimum_density - _density_tolerance)
    mooseError(name(),
               ": implicit temporary-density solve lost positivity. Minimum density is ",
               global_minimum,
               ", below the admissible minimum ",
               _minimum_density,
               ".");
  if (global_maximum > _maximum_density + _density_tolerance)
    mooseError(name(),
               ": implicit temporary-density solve lost boundedness. Maximum density is ",
               global_maximum,
               ", above the admissible maximum ",
               _maximum_density,
               ".");

  for (auto face = _fe_problem.mesh().ownedFaceInfoBegin();
       face != _fe_problem.mesh().ownedFaceInfoEnd();
       ++face)
    _mass_flux[(*face)->id()] = _mass_advection_kernel->faceMassFluxDensity(**face);
}

void
LowMachImplicitMassFlux::synchronizeDensityToEOS()
{
  if (!_mass_advection_kernel)
    mooseError(name(), ": mass-advection kernel has not been linked.");

  const auto & variable = _mass_advection_kernel->variable();
  auto & system = _fe_problem.getLinearSystem(variable.sys().number());
  auto & solution = system.solution();
  auto & current_local_solution = *system.system().current_local_solution;
  const auto system_number = system.number();
  const auto variable_number = variable.number();
  const auto state = Moose::currentState();

  for (const auto & elem_info : _fe_problem.mesh().elemInfoVector())
  {
    if (!_mass_advection_kernel->hasBlocks(elem_info->subdomain_id()))
      continue;

    const dof_id_type dof = elem_info->dofIndices()[system_number][variable_number];
    if (dof == DofObject::invalid_id || dof < solution.first_local_index() ||
        dof >= solution.last_local_index())
      continue;

    const Real eos_density = _eos_density(makeElemArg(elem_info->elem()), state);
    if (!std::isfinite(eos_density) || eos_density < _minimum_density - _density_tolerance ||
        eos_density > _maximum_density + _density_tolerance)
      mooseError(name(),
                 ": EOS synchronization produced inadmissible density ",
                 eos_density,
                 " at density DOF ",
                 dof,
                 ".");

    solution.set(dof, eos_density);
    current_local_solution.set(dof, eos_density);
  }

  solution.close();
  current_local_solution.close();
  system.setSolution(current_local_solution);
}
