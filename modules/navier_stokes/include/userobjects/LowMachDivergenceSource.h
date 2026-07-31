//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "BlockRestrictable.h"
#include "CellCenteredMapFunctor.h"
#include "GeneralUserObject.h"
#include "NonADFunctorInterface.h"

#include <unordered_map>
#include <vector>

class LinearFVLowMachEnthalpyAdvection;
class TimeIntegrator;

/**
 * Publishes the low-Mach velocity-divergence source from the converged enthalpy balance.
 *
 * The cell heating rate is reconstructed from the exact discrete conservative enthalpy LHS,
 * including the frozen face flux used by the final temperature equation. This avoids evaluating
 * thermal diffusion through a second, inconsistent stencil.
 */
class LowMachDivergenceSource : public GeneralUserObject,
                                public NonADFunctorInterface,
                                public BlockRestrictable
{
public:
  static InputParameters validParams();

  LowMachDivergenceSource(const InputParameters & params);

  void initialSetup() override;
  void meshChanged() override;
  void initialize() override {}
  void execute() override {}
  void finalize() override {}

  void beginFixedPointEnthalpySolve();
  void finishFixedPointEnthalpyInitialization();
  void beginEnthalpyAssembly();
  void publish();

  const SolverSystemName & systemName() const { return _system_name; }
  const MooseFunctorName & sourceName() const { return _source_name; }
  const MooseFunctorName & materialFractionName() const { return _material_fraction_name; }

protected:
  void linkEnthalpyAdvectionKernel();
  void initializeSourceStorage();

  const SolverSystemName _system_name;
  const std::string _enthalpy_advection_kernel_name;
  const Moose::Functor<Real> & _temporary_density;
  const Moose::Functor<Real> & _eos_density;
  const Moose::Functor<Real> & _specific_enthalpy;
  const Moose::Functor<Real> & _dliquid_fraction_dh;
  const MooseFunctorName _material_fraction_name;
  const Moose::Functor<Real> & _material_fraction;
  const Moose::Functor<Real> & _rho_solid;
  const Moose::Functor<Real> & _rho_liquid;
  const Moose::Functor<Real> & _h_solid;
  const Moose::Functor<Real> & _h_liquid;
  const TimeIntegrator * _time_integrator = nullptr;

  const MooseFunctorName _source_name;
  CellCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _source;
  std::vector<LinearFVLowMachEnthalpyAdvection *> _enthalpy_advection_kernels;
};
