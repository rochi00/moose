//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "LinearFVLowMachSolidDrag.h"

registerMooseObject("NavierStokesApp", LinearFVLowMachSolidDrag);

InputParameters
LinearFVLowMachSolidDrag::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Adds the Carman-Kozeny solid-phase penalization implicitly to a momentum component.");
  params.addRequiredParam<MooseFunctorName>("drag_coefficient",
                                            "Carman-Kozeny momentum drag coefficient.");
  return params;
}

LinearFVLowMachSolidDrag::LinearFVLowMachSolidDrag(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _coefficient_name(getParam<MooseFunctorName>("drag_coefficient")),
    _coefficient(getFunctor<Real>("drag_coefficient"))
{
}

Real
LinearFVLowMachSolidDrag::computeMatrixContribution()
{
  return _coefficient(makeElemArg(_current_elem_info->elem()), determineState()) *
         _current_elem_volume;
}

Real
LinearFVLowMachSolidDrag::computeRightHandSideContribution()
{
  return 0.0;
}
