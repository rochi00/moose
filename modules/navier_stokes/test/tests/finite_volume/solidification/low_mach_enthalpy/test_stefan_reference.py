import math
import unittest

from stefan_reference import CASES, interface_state, stefan_residual, temperature


class TestStefanReference(unittest.TestCase):
    def test_documented_interface_states_at_ten_seconds(self):
        expected = {
            "matched": (1.12533235, 0.0422695, 0.0),
            "expansion": (2.55867445, 0.0920167, 0.00374883),
            "shrinkage": (0.48781783, 0.0407668, -0.00896869),
        }
        for name, (lam, position, velocity) in expected.items():
            with self.subTest(case=name):
                state = interface_state(CASES[name], 10.0, include_kinetic_energy=True)
                self.assertAlmostEqual(state["lambda"], lam, places=7)
                self.assertAlmostEqual(state["position"], position, places=7)
                self.assertAlmostEqual(state["liquid_velocity"], velocity, places=7)

    def test_enthalpy_compatible_root_satisfies_energy_balance(self):
        for problem in CASES.values():
            state = interface_state(problem, 10.0)
            residual = stefan_residual(problem, state["lambda"], 10.0)
            heat_flux_scale = (
                problem.rho_solid
                * problem.effective_latent_heat
                * state["lambda"]
                * math.sqrt(problem.diffusivity_liquid)
            )
            self.assertLess(abs(residual / heat_flux_scale), 1.0e-12)

    def test_temperature_boundary_and_interface_values(self):
        for problem in CASES.values():
            state = interface_state(problem, 10.0)
            self.assertAlmostEqual(temperature(problem, 0.0, 10.0), problem.temperature_cold)
            self.assertAlmostEqual(
                temperature(problem, state["position"], 10.0),
                problem.temperature_melt,
            )
            self.assertAlmostEqual(
                temperature(problem, 1.0, 10.0),
                problem.temperature_initial,
                places=12,
            )

    def test_liquid_velocity_obeys_mass_jump(self):
        for problem in CASES.values():
            state = interface_state(problem, 10.0)
            expected = (1.0 - problem.density_ratio) * state["speed"]
            self.assertEqual(state["liquid_velocity"], expected)


if __name__ == "__main__":
    unittest.main()
