//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KokkosLinearFVElementalKernel.h"

class TimeIntegrator;

/**
 * Kokkos linear finite volume time derivative with a constant multiplier,
 *   factor * du/dt,
 * discretized with the implicit Euler or BDF2 time integrator of the variable. It is the Kokkos
 * analog of LinearFVTimeDerivative with a constant factor.
 */
class KokkosLinearFVTimeDerivative : public Moose::Kokkos::LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  KokkosLinearFVTimeDerivative(const InputParameters & parameters);

  template <typename Derived>
  KOKKOS_FUNCTION CellContribution computeCellContribution(const CellDatum & datum) const;

protected:
  virtual void dispatchContributions() override;

private:
  /// Constant multiplier on the time derivative
  const Real _factor;
  /// Time integrator of the variable
  const TimeIntegrator & _time_integrator;
  /// Old and older cell values of the variable
  const Moose::Kokkos::LinearFVCellValues _u_old;
  const Moose::Kokkos::LinearFVCellValues _u_older;

  /// Coefficients of the current time step, set before each dispatch: the matrix entry is
  /// factor * V * _matrix_coefficient and the right-hand side
  /// factor * V * (_old_coefficient * u_old + _older_coefficient * u_older)
  ///@{
  Real _matrix_coefficient = 0;
  Real _old_coefficient = 0;
  Real _older_coefficient = 0;
  ///@}
};

template <typename Derived>
KOKKOS_FUNCTION KokkosLinearFVTimeDerivative::CellContribution
KokkosLinearFVTimeDerivative::computeCellContribution(const CellDatum & datum) const
{
  const auto scale = _factor * datum.elemVolume();
  auto rhs = _old_coefficient * _u_old(kokkosSystems(), datum.elem());
  if (_older_coefficient)
    rhs += _older_coefficient * _u_older(kokkosSystems(), datum.elem());
  return {scale * _matrix_coefficient, scale * rhs};
}
