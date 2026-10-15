//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVGreenGaussGradient.h"

#include "ComputeLinearFVGreenGaussGradientFaceThread.h"
#include "ComputeLinearFVGreenGaussGradientVolumeThread.h"
#include "FEProblemBase.h"
#include "SystemBase.h"

#ifdef MOOSE_KOKKOS_ENABLED
#include "KokkosGreenGaussGradient.h"
#include "LinearSystem.h"
#endif

registerMooseObject("MooseApp", FVGreenGaussGradient);

InputParameters
FVGreenGaussGradient::validParams()
{
  InputParameters params = FVGradientMethod::validParams();
  params.addClassDescription("Green-Gauss cell-centered gradient method.");
  return params;
}

FVGreenGaussGradient::FVGreenGaussGradient(const InputParameters & params)
  : FVGradientMethod(params)
{
}

FVGreenGaussGradient::~FVGreenGaussGradient() = default;

void
FVGreenGaussGradient::computeGradientWithoutLimiter(
    SystemBase & system,
    GradientContainer & gradient,
    const std::unordered_set<unsigned int> & all_variable_numbers) const
{
#ifdef MOOSE_KOKKOS_ENABLED
  // The boundary conditions of a variable whose linear FV objects are all Kokkos objects are
  // visible to the Kokkos computation only, so its gradient is computed there
  const auto [kokkos_variable_numbers, variable_numbers] =
      splitKokkosVariables(system, all_variable_numbers);

  if (!kokkos_variable_numbers.empty())
  {
    if (!_kokkos_gradient)
      _kokkos_gradient = std::make_unique<Moose::Kokkos::GreenGaussGradient>();
    _kokkos_gradient->compute(system, gradient, kokkos_variable_numbers);
  }
  if (variable_numbers.empty())
    return;
#else
  const auto & variable_numbers = all_variable_numbers;
#endif

  auto & fe_problem = system.feProblem();

  PARALLEL_TRY
  {
    using FaceInfoRange = ComputeLinearFVGreenGaussGradientFaceThread::FaceInfoRange;
    FaceInfoRange face_info_range(fe_problem.mesh().ownedFaceInfoBegin(),
                                  fe_problem.mesh().ownedFaceInfoEnd());

    ComputeLinearFVGreenGaussGradientFaceThread gradient_face_thread(
        fe_problem, system, gradient, variable_numbers);
    Threads::parallel_reduce(face_info_range, gradient_face_thread, system.numThreads());
  }
  fe_problem.checkExceptionAndStopSolve();

  for (auto & vec : gradient)
    vec->close();

  PARALLEL_TRY
  {
    using ElemInfoRange = ComputeLinearFVGreenGaussGradientVolumeThread::ElemInfoRange;
    ElemInfoRange elem_info_range(fe_problem.mesh().ownedElemInfoBegin(),
                                  fe_problem.mesh().ownedElemInfoEnd());

    ComputeLinearFVGreenGaussGradientVolumeThread gradient_volume_thread(
        fe_problem, system, gradient, variable_numbers);
    Threads::parallel_reduce(elem_info_range, gradient_volume_thread, system.numThreads());
  }
  fe_problem.checkExceptionAndStopSolve();
}

#ifdef MOOSE_KOKKOS_ENABLED
void
FVGreenGaussGradient::applyLimiter(
    SystemBase & system,
    GradientContainer & gradient,
    const std::unordered_set<unsigned int> & all_variable_numbers) const
{
  // The gradients computed with Kokkos are limited on the device too
  const auto [kokkos_variable_numbers, variable_numbers] =
      splitKokkosVariables(system, all_variable_numbers);

  if (!kokkos_variable_numbers.empty())
    _kokkos_gradient->limit(system, gradient, kokkos_variable_numbers, limiterType());
  if (!variable_numbers.empty())
    FVGradientMethod::applyLimiter(system, gradient, variable_numbers);
}

std::pair<std::unordered_set<unsigned int>, std::unordered_set<unsigned int>>
FVGreenGaussGradient::splitKokkosVariables(
    const SystemBase & system, const std::unordered_set<unsigned int> & variable_numbers)
{
  std::pair<std::unordered_set<unsigned int>, std::unordered_set<unsigned int>> split;
  const auto * const linear_system = dynamic_cast<const LinearSystem *>(&system);
  for (const auto variable_number : variable_numbers)
    if (linear_system && linear_system->hasOnlyKokkosLinearFVObjects(variable_number))
      split.first.insert(variable_number);
    else
      split.second.insert(variable_number);
  return split;
}
#endif
