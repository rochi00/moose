//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVGradientComponentTestAux.h"

#include "MooseLinearVariableFV.h"
#include "MooseMesh.h"

registerMooseObject("MooseTestApp", LinearFVGradientComponentTestAux);

InputParameters
LinearFVGradientComponentTestAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<VariableName>("linear_variable",
                                        "The linear FV variable whose gradient is output.");
  params.addRequiredParam<unsigned int>("component", "The gradient component.");
  params.addClassDescription(
      "Test-only output of one component of the cell gradient of a linear FV variable.");
  return params;
}

LinearFVGradientComponentTestAux::LinearFVGradientComponentTestAux(
    const InputParameters & parameters)
  : AuxKernel(parameters),
    _linear_variable([this]() -> MooseLinearVariableFV<Real> &
                     {
                       auto * const variable = dynamic_cast<MooseLinearVariableFV<Real> *>(
                           &_subproblem.getVariable(0, getParam<VariableName>("linear_variable")));
                       if (!variable)
                         paramError("linear_variable", "Must be a linear FV variable.");
                       return *variable;
                     }()),
    _component(getParam<unsigned int>("component"))
{
  if (isNodal())
    mooseError("This aux kernel needs an elemental variable.");
  _linear_variable.requestCellGradients();
}

Real
LinearFVGradientComponentTestAux::computeValue()
{
  return _linear_variable.gradSlnComponent(_mesh.elemInfo(_current_elem->id()), _component);
}
