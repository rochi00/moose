//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"
#include "KokkosFunctionParser.h"

/**
 * A property of the Kokkos linear finite volume objects computed from a parsed expression of
 * linear finite volume variables, numbers, space (x, y, z) and time (t). It is the Kokkos analog of
 * ParsedFunctorMaterial, with the same parameters: the Kokkos objects evaluate the expression on
 * the device, on cells with the cell values of the variables and on faces with their face values,
 * as the host functors. It defines no host functor.
 */
class KokkosParsedFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  KokkosParsedFunctorMaterial(const InputParameters & parameters);

  /// The maximum number of functors in the expression
  static constexpr unsigned int max_functors = 8;

  /// Get the name of the property
  const std::string & propertyName() const { return _property_name; }
  /// Get the functors of the expression, in the order of the arguments of the evaluator
  const std::vector<std::string> & functorNames() const { return _functor_names; }
  /// Get the evaluator of the expression, whose arguments are the values of functorNames()
  const Moose::Kokkos::RPNEvaluator & evaluator() const { return _evaluator; }

private:
  /// The functors of the expression
  const std::vector<std::string> _functor_names;
  /// The name of the property
  const std::string _property_name;
  /// The constants pi and e, as ParsedFunctorMaterial defines them
  ///@{
  const Real _pi;
  const Real _e;
  ///@}
  /// The evaluator of the expression
  Moose::Kokkos::RPNEvaluator _evaluator;
};
