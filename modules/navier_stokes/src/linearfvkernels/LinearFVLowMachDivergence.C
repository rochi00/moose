//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVLowMachDivergence.h"

registerMooseObject("NavierStokesApp", LinearFVLowMachDivergence);

InputParameters
LinearFVLowMachDivergence::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Adds the discretely consistent phase-change velocity-divergence source to the low-Mach "
      "pressure equation.");
  params.addRequiredParam<MooseFunctorName>("divergence_source",
                                            "Cell low-Mach divergence-source functor.");
  return params;
}

LinearFVLowMachDivergence::LinearFVLowMachDivergence(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _source_name(getParam<MooseFunctorName>("divergence_source")),
    _source(getFunctor<Real>("divergence_source"))
{
}

Real
LinearFVLowMachDivergence::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVLowMachDivergence::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  return _source(elem_arg, determineState()) * _current_elem_volume;
}
