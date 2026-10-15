//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVGradientMethod.h"

#include <map>
#include <memory>

namespace Moose::Kokkos
{

/**
 * Green-Gauss cell-centered gradients of linear finite-volume variables computed on device, as
 * FVGreenGaussGradient: face values come from the linear interpolation of the two cell values,
 * from the value relation of the Kokkos boundary condition on boundary faces, or from the cell
 * value on boundary faces without one. Each cell sums its sides in a fixed order, so the result is
 * identical in every run. Used by KokkosFVGreenGaussGradient, and by FVGreenGaussGradient for the
 * variables whose linear finite volume objects are all Kokkos objects.
 */
class GreenGaussGradient
{
public:
  GreenGaussGradient();
  ~GreenGaussGradient();

  /**
   * Compute the gradients without limiter
   * @param system The system of the variables
   * @param gradient The component vectors to write the gradients into
   * @param variable_numbers The variables
   */
  void compute(SystemBase & system,
               FVGradientMethod::GradientContainer & gradient,
               const std::unordered_set<unsigned int> & variable_numbers);

  /**
   * Limit gradients computed by compute(), as FVGradientMethod::applyLimiter()
   * @param system The system of the variables
   * @param gradient The component vectors holding the gradients, limited in place
   * @param variable_numbers The variables
   * @param limiter_type The limiter
   */
  void limit(SystemBase & system,
             FVGradientMethod::GradientContainer & gradient,
             const std::unordered_set<unsigned int> & variable_numbers,
             Moose::FV::GradientLimiterType limiter_type);

private:
  /// Device tables of one variable, built at its first gradient computation
  struct VariableData;

  /**
   * Get the device tables of a variable, building them if needed
   * @param system The variable's system
   * @param variable_number The variable number in its system
   */
  VariableData & variableData(SystemBase & system, unsigned int variable_number);

  /// Device tables by system and variable number
  std::map<std::pair<unsigned int, unsigned int>, std::unique_ptr<VariableData>> _variable_data;
};

} // namespace Moose::Kokkos
