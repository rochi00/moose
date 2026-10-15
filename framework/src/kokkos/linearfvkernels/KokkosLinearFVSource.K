//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosLinearFVSource.h"

registerKokkosLinearFVKernel("MooseApp", KokkosLinearFVSource);

InputParameters
KokkosLinearFVSource::validParams()
{
  InputParameters params = Moose::Kokkos::LinearFVElementalKernel::validParams();
  params.addParam<MooseFunctorName>("source_density", 1.0, "The source density.");
  params.addParam<MooseFunctorName>(
      "scaling_factor", 1.0, "Scaling factor for the source term (for example a density).");
  params.addClassDescription("Kokkos LinearFV source kernel.");
  return params;
}

KokkosLinearFVSource::KokkosLinearFVSource(const InputParameters & parameters)
  : Moose::Kokkos::LinearFVElementalKernel(parameters)
{
  PropertyUse use;
  use.cells = true;
  requestLinearFVProperty("source_density", _source_density, use);
  requestLinearFVProperty("scaling_factor", _scale, use);
}
