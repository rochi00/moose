//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#pragma once

#include "BlockRestrictable.h"
#include "GeneralUserObject.h"
#include "MooseLinearVariableFV.h"
#include "NonADFunctorInterface.h"

#include <unordered_map>

class SystemBase;

/**
 * Restores exact thermodynamic consistency after a Newton-linearized low-Mach enthalpy solve.
 *
 * Before the linear temperature solve, beginNewtonIteration() stores T^m, h^m, dh/dT^m, and
 * liquid fraction. After the solve, restoreExactState() first forms
 *
 *   h^(m+1) = h^m + (dh/dT)^m (T_linear^(m+1) - T^m)
 *
 * and then replaces the linear temperature with the exact inverse T(h^(m+1)). A damped update
 * interpolates in enthalpy, not temperature, so every line-search trial remains on the exact
 * constitutive curve. The returned correction is the global discrete L2 liquid-fraction update
 * specified by the low-Mach enthalpy algorithm.
 */
class LowMachEnthalpyNewtonState : public GeneralUserObject,
                                   public NonADFunctorInterface,
                                   public BlockRestrictable
{
public:
  static InputParameters validParams();

  LowMachEnthalpyNewtonState(const InputParameters & params);

  void initialize() override {}
  void execute() override {}
  void finalize() override {}

  void beginNewtonIteration();
  void captureNewtonUpdate();
  Real restoreExactState(Real step_length = 1.0);

  const SolverSystemName & systemName() const { return _system_name; }
  const MooseFunctorName & materialFractionName() const { return _material_fraction_name; }

protected:
  struct IterationState
  {
    Real temperature;
    Real enthalpy;
    Real dh_dT;
    Real liquid_fraction;
    Real material_fraction;
    Real cp_gas;
    Real cp_solid;
    Real cp_liquid;
    Real rho_solid;
    Real rho_liquid;
    Real h_solid;
    Real h_liquid;
    Real T_solidus;
    Real T_liquidus;
    Real T_reference;
    Real full_step_enthalpy;
  };

  Real exactTemperature(Real enthalpy, const IterationState & state) const;
  Real exactLiquidFraction(Real enthalpy, const IterationState & state) const;

  MooseLinearVariableFVReal & _temperature_variable;
  SystemBase & _temperature_system;
  const SolverSystemName _system_name;
  const unsigned int _system_number;
  const unsigned int _variable_number;

  const MooseFunctorName _material_fraction_name;
  const Moose::Functor<Real> & _material_fraction;
  const Moose::Functor<Real> & _specific_enthalpy;
  const Moose::Functor<Real> & _dh_dT;
  const Moose::Functor<Real> & _liquid_fraction;
  const Moose::Functor<Real> & _cp_gas;
  const Moose::Functor<Real> & _cp_solid;
  const Moose::Functor<Real> & _cp_liquid;
  const Moose::Functor<Real> & _rho_solid;
  const Moose::Functor<Real> & _rho_liquid;
  const Moose::Functor<Real> & _h_solid;
  const Moose::Functor<Real> & _h_liquid;
  const Moose::Functor<Real> & _T_solidus;
  const Moose::Functor<Real> & _T_liquidus;
  const Moose::Functor<Real> & _T_reference;

  std::unordered_map<dof_id_type, IterationState> _iteration_state;
  bool _has_captured_update = false;
};
