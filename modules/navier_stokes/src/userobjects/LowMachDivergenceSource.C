//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LowMachDivergenceSource.h"

#include "FVUtils.h"
#include "LinearFVLowMachEnthalpyAdvection.h"
#include "LowMachEnthalpyThermodynamics.h"
#include "TheWarehouse.h"
#include "TimeIntegrator.h"

#include "timpi/parallel_sync.h"

#include <cmath>
#include <map>

registerMooseObject("NavierStokesApp", LowMachDivergenceSource);

InputParameters
LowMachDivergenceSource::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params += BlockRestrictable::validParams();
  params.addClassDescription(
      "Publishes the low-Mach divergence source using the exact converged conservative enthalpy "
      "balance and its assembled common-mass-flux face operator.");
  params.addRequiredParam<SolverSystemName>("system", "Low-Mach temperature system.");
  params.addRequiredParam<std::string>("enthalpy_advection_kernel",
                                       "Name of the LinearFVLowMachEnthalpyAdvection kernel.");
  params.addRequiredParam<MooseFunctorName>(
      "temporary_density", "End-of-step density produced by the conservative mass solve.");
  params.addRequiredParam<MooseFunctorName>("eos_density", "Current thermodynamic EOS density.");
  params.addRequiredParam<MooseFunctorName>("specific_enthalpy",
                                            "Current exact specific enthalpy.");
  params.addRequiredParam<MooseFunctorName>(
      "dliquid_fraction_dh",
      "Exact current derivative of liquid fraction with respect to enthalpy.");
  params.addRequiredParam<MooseFunctorName>("material_fraction", "Gas/PCM material fraction.");
  params.addRequiredParam<MooseFunctorName>("rho_solid", "Solid PCM density.");
  params.addRequiredParam<MooseFunctorName>("rho_liquid", "Liquid PCM density.");
  params.addRequiredParam<MooseFunctorName>("h_solid", "Solidus specific enthalpy.");
  params.addRequiredParam<MooseFunctorName>("h_liquid", "Liquidus specific enthalpy.");
  params.addParam<MooseFunctorName>(
      "source_name", "low_mach_divergence_source", "Published cell divergence-source functor.");

  ExecFlagEnum & exec_enum = params.set<ExecFlagEnum>("execute_on", true);
  exec_enum.addAvailableFlags(EXEC_NONE);
  exec_enum = {EXEC_NONE};
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

LowMachDivergenceSource::LowMachDivergenceSource(const InputParameters & params)
  : GeneralUserObject(params),
    NonADFunctorInterface(this),
    BlockRestrictable(this),
    _system_name(getParam<SolverSystemName>("system")),
    _enthalpy_advection_kernel_name(getParam<std::string>("enthalpy_advection_kernel")),
    _temporary_density(getFunctor<Real>("temporary_density")),
    _eos_density(getFunctor<Real>("eos_density")),
    _specific_enthalpy(getFunctor<Real>("specific_enthalpy")),
    _dliquid_fraction_dh(getFunctor<Real>("dliquid_fraction_dh")),
    _material_fraction_name(getParam<MooseFunctorName>("material_fraction")),
    _material_fraction(getFunctor<Real>("material_fraction")),
    _rho_solid(getFunctor<Real>("rho_solid")),
    _rho_liquid(getFunctor<Real>("rho_liquid")),
    _h_solid(getFunctor<Real>("h_solid")),
    _h_liquid(getFunctor<Real>("h_liquid")),
    _source_name(getParam<MooseFunctorName>("source_name")),
    _source(_fe_problem.mesh(), blockIDs(), _source_name, true)
{
  for (const auto tid : make_range(libMesh::n_threads()))
    UserObject::_subproblem.addFunctor(_source_name, _source, tid);
}

void
LowMachDivergenceSource::initialSetup()
{
  linkEnthalpyAdvectionKernel();
  initializeSourceStorage();
}

void
LowMachDivergenceSource::meshChanged()
{
  linkEnthalpyAdvectionKernel();
  initializeSourceStorage();
}

void
LowMachDivergenceSource::linkEnthalpyAdvectionKernel()
{
  const auto system_number = _fe_problem.linearSysNum(_system_name);
  std::vector<LinearFVFluxKernel *> kernels;
  _fe_problem.theWarehouse()
      .query()
      .condition<AttribSysNum>(system_number)
      .condition<AttribSystem>("LinearFVFluxKernel")
      .condition<AttribName>(_enthalpy_advection_kernel_name)
      .queryInto(kernels);

  if (kernels.empty())
    paramError("enthalpy_advection_kernel",
               "No linear FV flux kernel named '",
               _enthalpy_advection_kernel_name,
               "' in system '",
               _system_name,
               "' was found.");

  _enthalpy_advection_kernels.clear();
  _time_integrator = nullptr;
  for (auto * kernel : kernels)
  {
    auto * enthalpy_advection = dynamic_cast<LinearFVLowMachEnthalpyAdvection *>(kernel);
    if (!enthalpy_advection)
      paramError("enthalpy_advection_kernel",
                 "Every thread-local kernel named '",
                 _enthalpy_advection_kernel_name,
                 "' must be a LinearFVLowMachEnthalpyAdvection.");

    const auto & variable = enthalpy_advection->variable();
    const auto * const time_integrator = &variable.sys().getTimeIntegrator(variable.number());
    if (_time_integrator && _time_integrator != time_integrator)
      paramError("enthalpy_advection_kernel",
                 "All thread-local enthalpy-advection kernels must use the same time integrator.");
    _time_integrator = time_integrator;
    _enthalpy_advection_kernels.push_back(enthalpy_advection);
  }
}

void
LowMachDivergenceSource::initializeSourceStorage()
{
  _source.clear();
  for (const auto * elem_info : _fe_problem.mesh().elemInfoVector())
    _source[elem_info->elem()->id()] = 0.0;
}

void
LowMachDivergenceSource::beginFixedPointEnthalpySolve()
{
  if (_enthalpy_advection_kernels.empty())
    mooseError(name(), ": enthalpy-advection kernels have not been linked.");
  for (auto * kernel : _enthalpy_advection_kernels)
    kernel->beginFixedPointFluxCapture();
}

void
LowMachDivergenceSource::finishFixedPointEnthalpyInitialization()
{
  for (auto * kernel : _enthalpy_advection_kernels)
    kernel->finishFixedPointFluxCapture();
}

void
LowMachDivergenceSource::beginEnthalpyAssembly()
{
  if (_enthalpy_advection_kernels.empty())
    mooseError(name(), ": enthalpy-advection kernels have not been linked.");
  for (auto * kernel : _enthalpy_advection_kernels)
    kernel->clearAssembledFaceEnthalpyFluxes();
}

void
LowMachDivergenceSource::publish()
{
  if (_enthalpy_advection_kernels.empty())
    mooseError(name(), ": enthalpy-advection kernels have not been linked.");
  if (!_time_integrator)
    mooseError(name(), ": no enthalpy time integrator has been linked.");
  if (_dt <= 0.0)
    mooseError(name(), ": a positive time step is required to construct the divergence source.");

  using FaceKey = std::pair<dof_id_type, unsigned int>;
  const auto face_key = [](const FaceInfo & face_info)
  { return FaceKey(face_info.elem().id(), face_info.elemSideID()); };

  std::map<FaceKey, Real> assembled_face_enthalpy_flux;
  for (const auto * face_info : _fe_problem.mesh().faceInfo())
    for (const auto * kernel : _enthalpy_advection_kernels)
      if (kernel->hasAssembledFaceEnthalpyFlux(*face_info))
      {
        assembled_face_enthalpy_flux.emplace(face_key(*face_info),
                                             kernel->assembledFaceEnthalpyFlux(*face_info));
        break;
      }

  if (n_processors() > 1)
  {
    using Datum = std::pair<FaceKey, Real>;
    std::unordered_map<processor_id_type, std::vector<Datum>> push_data;
    for (const auto * face_info : _fe_problem.mesh().faceInfo())
    {
      const auto * neighbor = face_info->neighborPtr();
      if (face_info->processor_id() != processor_id() || !neighbor ||
          neighbor->processor_id() == processor_id())
        continue;

      const auto key = face_key(*face_info);
      const auto flux = assembled_face_enthalpy_flux.find(key);
      if (flux == assembled_face_enthalpy_flux.end())
        mooseError(name(),
                   ": the owning rank did not assemble enthalpy flux on partition face ",
                   face_info->id(),
                   " in the final Newton equation.");
      push_data[neighbor->processor_id()].emplace_back(key, flux->second);
    }

    auto receive_partition_fluxes =
        [&assembled_face_enthalpy_flux](const processor_id_type,
                                        const std::vector<Datum> & received_fluxes)
    {
      for (const auto & [key, flux] : received_fluxes)
        assembled_face_enthalpy_flux.emplace(key, flux);
    };
    TIMPI::push_parallel_vector_data(_communicator, push_data, receive_partition_fluxes);
  }

  const auto state = Moose::currentState();
  const auto coefficients = _time_integrator->timeDerivativeCoefficients();
  const bool backward_euler =
      coefficients.size() == 2 && coefficients[0] == 1.0 && coefficients[1] == -1.0;
  for (const auto * elem_info : _fe_problem.mesh().elemInfoVector())
  {
    const Elem * const elem = elem_info->elem();
    if (!hasBlocks(elem->subdomain_id()) || elem->processor_id() != processor_id())
      continue;

    const auto elem_arg = makeElemArg(elem);
    Real net_outward_enthalpy_flux = 0.0;
    for (const auto side : make_range(elem->n_sides()))
    {
      const Elem * const neighbor = elem->neighbor_ptr(side);
      const bool elem_owns_face = Moose::FV::elemHasFaceInfo(*elem, neighbor);
      const FaceInfo * const face_info =
          elem_owns_face
              ? _fe_problem.mesh().faceInfo(elem, side)
              : _fe_problem.mesh().faceInfo(neighbor, neighbor->which_neighbor_am_i(elem));
      if (!face_info)
        mooseError(
            name(), ": no FaceInfo is available for element ", elem->id(), " side ", side, ".");

      const Real orientation = face_info->elemPtr() == elem ? 1.0 : -1.0;
      const auto assembled_flux = assembled_face_enthalpy_flux.find(face_key(*face_info));
      if (assembled_flux == assembled_face_enthalpy_flux.end())
        mooseError(name(),
                   ": no local or synchronized enthalpy flux is available for face ",
                   face_info->id(),
                   " in the final Newton equation.");
      net_outward_enthalpy_flux += orientation * assembled_flux->second;
    }

    const Real cell_volume = elem_info->volume() * elem_info->coordFactor();
    const Real current_enthalpy = _specific_enthalpy(elem_arg, state);
    Real enthalpy_time_derivative =
        coefficients.front() * _temporary_density(elem_arg, state) * current_enthalpy;
    for (const auto state_index : index_range(coefficients))
    {
      if (state_index == 0)
        continue;

      const Moose::StateArg history_state(state_index, Moose::SolutionIterationType::Time);
      enthalpy_time_derivative += coefficients[state_index] *
                                  _eos_density(elem_arg, history_state) *
                                  _specific_enthalpy(elem_arg, history_state);
    }
    enthalpy_time_derivative /= _dt;
    const Real heating_rate = enthalpy_time_derivative + net_outward_enthalpy_flux / cell_volume;

    const Real eos_density = _eos_density(elem_arg, state);
    const Real material_fraction = _material_fraction(elem_arg, state);
    const Real rho_solid = _rho_solid(elem_arg, state);
    const Real rho_liquid = _rho_liquid(elem_arg, state);
    const Real source =
        backward_euler ? LowMachEnthalpy::backwardEulerDivergenceSource(material_fraction,
                                                                        eos_density,
                                                                        rho_solid,
                                                                        rho_liquid,
                                                                        current_enthalpy,
                                                                        _h_solid(elem_arg, state),
                                                                        _h_liquid(elem_arg, state),
                                                                        heating_rate,
                                                                        _dt)
                       : LowMachEnthalpy::divergenceSource(material_fraction,
                                                           eos_density,
                                                           rho_solid,
                                                           rho_liquid,
                                                           _dliquid_fraction_dh(elem_arg, state),
                                                           heating_rate);

    if (!std::isfinite(source))
      mooseError(name(), ": non-finite divergence source produced in element ", elem->id(), ".");
    _source[elem->id()] = source;
  }
}
