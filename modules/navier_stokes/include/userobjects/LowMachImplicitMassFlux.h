//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "FaceCenteredMapFunctor.h"
#include "GeneralUserObject.h"
#include "NonADFunctorInterface.h"

#include <unordered_map>

class LinearFVLowMachMassAdvection;

/**
 * Publishes the converged face flux from LinearFVLowMachMassAdvection.
 *
 * The owning advection kernel supplies the reconstruction so the published flux is algebraically
 * identical to the one used by the conservative temporary-density equation.
 */
class LowMachImplicitMassFlux : public GeneralUserObject, public NonADFunctorInterface
{
public:
  static InputParameters validParams();

  LowMachImplicitMassFlux(const InputParameters & params);

  void initialSetup() override;
  void meshChanged() override;
  void initialize() override {}
  void execute() override;
  void finalize() override {}

  void publishConvergedMassFlux();
  void synchronizeDensityToEOS();

  const SolverSystemName & systemName() const { return _system_name; }
  const MooseFunctorName & massFluxName() const { return _mass_flux_name; }
  const MooseFunctorName & densityVariableName() const { return _density_variable_name; }

protected:
  void linkMassAdvectionKernel();
  void initializeFluxStorage();

  const SolverSystemName _system_name;
  const std::string _mass_advection_kernel_name;
  const MooseFunctorName _density_variable_name;
  const Real _minimum_density;
  const Real _maximum_density;
  const Real _density_tolerance;
  const Moose::Functor<Real> & _eos_density;

  const MooseFunctorName _mass_flux_name;
  FaceCenteredMapFunctor<Real, std::unordered_map<dof_id_type, Real>> _mass_flux;
  LinearFVLowMachMassAdvection * _mass_advection_kernel = nullptr;
};
