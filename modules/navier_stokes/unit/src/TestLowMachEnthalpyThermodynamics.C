//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "LowMachCUI.h"
#include "LowMachEnthalpyThermodynamics.h"

namespace
{
const LowMachEnthalpy::Parameters parameters{2.0, 3.0, 5.0, 0.1, 8.0, 6.0, 100.0, 10.0, 14.0, 0.0};
}

TEST(LowMachEnthalpyThermodynamicsTest, phaseEnthalpyBounds)
{
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::solidEnthalpy(parameters), 30.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::liquidEnthalpy(parameters), 146.0);
}

TEST(LowMachEnthalpyThermodynamicsTest, gasBranch)
{
  constexpr Real material_fraction = 0.49;
  constexpr Real temperature = 25.0;
  const Real enthalpy =
      LowMachEnthalpy::enthalpyFromTemperature(temperature, material_fraction, parameters);

  EXPECT_DOUBLE_EQ(enthalpy, 50.0);
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::temperatureFromEnthalpy(enthalpy, material_fraction, parameters),
      temperature);
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::liquidFractionFromEnthalpy(enthalpy, material_fraction, parameters), 0.0);
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::dEnthalpyDTemperature(temperature, material_fraction, parameters), 2.0);
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::dLiquidFractionDEnthalpy(enthalpy, material_fraction, parameters), 0.0);
}

TEST(LowMachEnthalpyThermodynamicsTest, pcmBranches)
{
  constexpr Real material_fraction = 1.0;

  EXPECT_DOUBLE_EQ(LowMachEnthalpy::enthalpyFromTemperature(5.0, material_fraction, parameters),
                   15.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::temperatureFromEnthalpy(15.0, material_fraction, parameters),
                   5.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::liquidFractionFromEnthalpy(15.0, material_fraction, parameters),
                   0.0);

  constexpr Real mushy_enthalpy = 88.0;
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::temperatureFromEnthalpy(mushy_enthalpy, material_fraction, parameters),
      12.0);
  EXPECT_NEAR(
      LowMachEnthalpy::liquidFractionFromEnthalpy(mushy_enthalpy, material_fraction, parameters),
      4.0 / 7.0,
      1e-14);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::dEnthalpyDTemperature(12.0, material_fraction, parameters),
                   29.0);

  EXPECT_DOUBLE_EQ(LowMachEnthalpy::temperatureFromEnthalpy(166.0, material_fraction, parameters),
                   18.0);
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::liquidFractionFromEnthalpy(166.0, material_fraction, parameters), 1.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::dEnthalpyDTemperature(18.0, material_fraction, parameters),
                   5.0);
}

TEST(LowMachEnthalpyThermodynamicsTest, exactInverse)
{
  for (const Real material_fraction : {0.0, 1.0})
    for (const Real temperature : {-2.0, 5.0, 10.0, 11.0, 12.0, 14.0, 18.0})
    {
      const Real enthalpy =
          LowMachEnthalpy::enthalpyFromTemperature(temperature, material_fraction, parameters);
      EXPECT_NEAR(LowMachEnthalpy::temperatureFromEnthalpy(enthalpy, material_fraction, parameters),
                  temperature,
                  1e-14);
    }
}

TEST(LowMachEnthalpyThermodynamicsTest, newtonUpdateRestoresExactMushyState)
{
  constexpr Real material_fraction = 1.0;
  constexpr Real old_temperature = 5.0;
  const Real old_enthalpy =
      LowMachEnthalpy::enthalpyFromTemperature(old_temperature, material_fraction, parameters);
  const Real old_dh_dT =
      LowMachEnthalpy::dEnthalpyDTemperature(old_temperature, material_fraction, parameters);
  const Real updated_enthalpy = LowMachEnthalpy::newtonUpdatedEnthalpy(
      old_enthalpy, old_dh_dT, old_temperature, 29.333333333333);

  EXPECT_NEAR(updated_enthalpy, 88.0, 1e-12);
  EXPECT_NEAR(
      LowMachEnthalpy::temperatureFromEnthalpy(updated_enthalpy, material_fraction, parameters),
      12.0,
      1e-12);
  EXPECT_NEAR(
      LowMachEnthalpy::liquidFractionFromEnthalpy(updated_enthalpy, material_fraction, parameters),
      4.0 / 7.0,
      1e-12);
}

TEST(LowMachEnthalpyThermodynamicsTest, liquidFractionDerivative)
{
  constexpr Real enthalpy = 88.0;
  constexpr Real material_fraction = 1.0;
  constexpr Real epsilon = 1e-6;
  const Real finite_difference = (LowMachEnthalpy::liquidFractionFromEnthalpy(
                                      enthalpy + epsilon, material_fraction, parameters) -
                                  LowMachEnthalpy::liquidFractionFromEnthalpy(
                                      enthalpy - epsilon, material_fraction, parameters)) /
                                 (2.0 * epsilon);

  EXPECT_NEAR(LowMachEnthalpy::dLiquidFractionDEnthalpy(enthalpy, material_fraction, parameters),
              finite_difference,
              1e-10);
}

TEST(LowMachEnthalpyThermodynamicsTest, eosAndLatentEnergy)
{
  constexpr Real material_fraction = 0.75;
  constexpr Real liquid_fraction = 4.0 / 7.0;

  EXPECT_NEAR(LowMachEnthalpy::density(material_fraction, liquid_fraction, parameters),
              5.167857142857143,
              1e-14);
  EXPECT_NEAR(
      LowMachEnthalpy::volumetricLatentEnthalpy(material_fraction, liquid_fraction, parameters),
      257.14285714285717,
      1e-13);
}

TEST(LowMachEnthalpyThermodynamicsTest, synchronizedThreePhaseProperty)
{
  EXPECT_NEAR(LowMachEnthalpy::threePhaseProperty(0.75, 4.0 / 7.0, 2.0, 3.0, 5.0),
              3.607142857142857,
              1e-14);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::threePhaseProperty(0.0, 1.0, 2.0, 3.0, 5.0), 2.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::threePhaseProperty(1.0, 0.0, 2.0, 3.0, 5.0), 3.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::threePhaseProperty(1.0, 1.0, 2.0, 3.0, 5.0), 5.0);
}

TEST(LowMachEnthalpyThermodynamicsTest, divergenceSource)
{
  const Real derivative = LowMachEnthalpy::dLiquidFractionDEnthalpy(88.0, 1.0, parameters);
  const Real eos_density = LowMachEnthalpy::density(1.0, 4.0 / 7.0, parameters);
  const Real expected = 2.0 * derivative * 584.0 / (eos_density * eos_density);

  EXPECT_NEAR(LowMachEnthalpy::divergenceSource(1.0, eos_density, 8.0, 6.0, derivative, 584.0),
              expected,
              1e-14);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::divergenceSource(0.0, 0.1, 8.0, 6.0, derivative, 584.0), 0.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::divergenceSource(1.0, 8.0, 8.0, 8.0, derivative, 584.0), 0.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::divergenceSource(1.0, 8.0, 8.0, 6.0, 0.0, 584.0), 0.0);
}

TEST(LowMachEnthalpyThermodynamicsTest, backwardEulerDivergenceMatchesFiniteEOSChange)
{
  constexpr Real material_fraction = 0.75;
  constexpr Real previous_enthalpy = 88.0;
  constexpr Real current_enthalpy = 15.0;
  constexpr Real dt = 0.25;
  const Real previous_liquid_fraction =
      LowMachEnthalpy::liquidFractionFromEnthalpy(previous_enthalpy, material_fraction, parameters);
  const Real current_liquid_fraction =
      LowMachEnthalpy::liquidFractionFromEnthalpy(current_enthalpy, material_fraction, parameters);
  const Real previous_density =
      LowMachEnthalpy::density(material_fraction, previous_liquid_fraction, parameters);
  const Real current_density =
      LowMachEnthalpy::density(material_fraction, current_liquid_fraction, parameters);
  const Real heating_rate = current_density * (current_enthalpy - previous_enthalpy) / dt;

  const Real source =
      LowMachEnthalpy::backwardEulerDivergenceSource(material_fraction,
                                                     current_density,
                                                     parameters.rho_solid,
                                                     parameters.rho_liquid,
                                                     current_enthalpy,
                                                     LowMachEnthalpy::solidEnthalpy(parameters),
                                                     LowMachEnthalpy::liquidEnthalpy(parameters),
                                                     heating_rate,
                                                     dt);

  EXPECT_NEAR(source, (previous_density / current_density - 1.0) / dt, 1e-14);
  EXPECT_NEAR(current_density * (1.0 + dt * source), previous_density, 1e-14);
}

TEST(LowMachEnthalpyThermodynamicsTest, backwardEulerDivergenceCapturesEndpointCrossing)
{
  constexpr Real material_fraction = 1.0;
  constexpr Real previous_enthalpy = 88.0;
  constexpr Real current_enthalpy = 15.0;
  constexpr Real dt = 0.25;
  const Real current_density = parameters.rho_solid;
  const Real heating_rate = current_density * (current_enthalpy - previous_enthalpy) / dt;

  EXPECT_LT(
      LowMachEnthalpy::backwardEulerDivergenceSource(material_fraction,
                                                     current_density,
                                                     parameters.rho_solid,
                                                     parameters.rho_liquid,
                                                     current_enthalpy,
                                                     LowMachEnthalpy::solidEnthalpy(parameters),
                                                     LowMachEnthalpy::liquidEnthalpy(parameters),
                                                     heating_rate,
                                                     dt),
      0.0);
  EXPECT_DOUBLE_EQ(
      LowMachEnthalpy::dLiquidFractionDEnthalpy(current_enthalpy, material_fraction, parameters),
      0.0);
}

TEST(LowMachEnthalpyThermodynamicsTest, solidDragCoefficient)
{
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::solidDragCoefficient(0.0, 0.0, 8.0, 1.0, 1e-3), 0.0);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::solidDragCoefficient(1.0, 1.0, 8.0, 1.0, 1e-3), 0.0);
  EXPECT_NEAR(LowMachEnthalpy::solidDragCoefficient(1.0, 4.0 / 7.0, 8.0, 1.0, 1e-3),
              504000.0 / 64343.0,
              1e-14);
  EXPECT_DOUBLE_EQ(LowMachEnthalpy::solidDragCoefficient(1.0, 0.0, 8.0, 1.0, 1e-3), 8000.0);
}

TEST(LowMachCUIReconstructionTest, korenCUIBranches)
{
  EXPECT_DOUBLE_EQ(LowMachCUI::faceValue(0.1, 0.0, 1.0), 0.3);
  EXPECT_DOUBLE_EQ(LowMachCUI::faceValue(0.5, 0.0, 1.0), 0.75);
  EXPECT_DOUBLE_EQ(LowMachCUI::faceValue(0.9, 0.0, 1.0), 1.0);
}

TEST(LowMachCUIReconstructionTest, korenCUIUsesDonorOutsideMonotoneRange)
{
  EXPECT_DOUBLE_EQ(LowMachCUI::faceValue(-0.2, 0.0, 1.0), -0.2);
  EXPECT_DOUBLE_EQ(LowMachCUI::faceValue(1.2, 0.0, 1.0), 1.2);
  EXPECT_DOUBLE_EQ(LowMachCUI::faceValue(4.0, 4.0, 4.0), 4.0);
}

TEST(LowMachCUIReconstructionTest, cubicCUIIsThirdOrderForSmoothCellAverages)
{
  const auto cell_average = [](const Real center, const Real h)
  { return std::exp(center) * 2.0 * std::sinh(h / 2.0) / h; };
  const auto face_error = [&cell_average](const Real h)
  {
    const Real far_upwind = cell_average(-h, h);
    const Real upwind = cell_average(0.0, h);
    const Real downwind = cell_average(h, h);
    return std::abs(LowMachCUI::faceValue(upwind, far_upwind, downwind) - std::exp(h / 2.0));
  };

  const Real coarse_error = face_error(0.2);
  const Real fine_error = face_error(0.1);
  const Real observed_order = std::log(coarse_error / fine_error) / std::log(2.0);

  EXPECT_GT(observed_order, 3.0);
  EXPECT_LT(observed_order, 3.1);
}
