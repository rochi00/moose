//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "KokkosParsedFunctorMaterial.h"

#include <cmath>

registerMooseObject("MooseApp", KokkosParsedFunctorMaterial);

InputParameters
KokkosParsedFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Computes a property of Kokkos linear finite volume objects from a "
                             "parsed expression of linear finite volume variables and numbers.");
  params.addRequiredCustomTypeParam<std::string>(
      "expression", "FunctionExpression", "Expression to parse for the new functor material");
  params.addParam<std::vector<std::string>>(
      "functor_names",
      {},
      "Functors to use in the parsed expression: linear finite volume variables or numbers");
  params.addParam<std::vector<std::string>>(
      "functor_symbols",
      {},
      "Symbolic name to use for each functor in 'functor_names' in the parsed expression. If not "
      "provided, then the actual functor names must be used in the parsed expression.");
  params.addRequiredParam<std::string>("property_name",
                                       "Name to give the new functor material property");
  return params;
}

KokkosParsedFunctorMaterial::KokkosParsedFunctorMaterial(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _functor_names(getParam<std::vector<std::string>>("functor_names")),
    _property_name(getParam<std::string>("property_name")),
    // As ParsedFunctorMaterial, which passes the constants as strings with six decimals
    _pi(std::stod(std::to_string(libMesh::pi))),
    _e(std::stod(std::to_string(std::exp(Real(1)))))
{
  if (_functor_names.size() > max_functors)
    paramError("functor_names", "At most ", max_functors, " functors are supported.");

  auto symbols = getParam<std::vector<std::string>>("functor_symbols");
  if (symbols.empty())
    symbols = _functor_names;
  else if (symbols.size() != _functor_names.size())
    paramError("functor_symbols",
               "The number of entries must be equal to either zero or the number of entries in "
               "'functor_names'.");

  Moose::Kokkos::RPNBuilder builder(getParam<std::string>("expression"), &_console);
  builder.addDefaultVariables();
  for (const auto & symbol : symbols)
    builder.addVariable(symbol);
  builder.addVariable("pi");
  builder.addVariable("e");
  builder.build();

  for (const auto i : index_range(symbols))
    builder.associateArgument(symbols[i], i);
  builder.associateScalar("pi", &_pi);
  builder.associateScalar("e", &_e);
  builder.finalize();

  _evaluator.init(builder);
}
