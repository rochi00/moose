import unittest

from stefan_reference import CASES, interface_state, temperature
from validate_stefan import (
    bulk_liquid_velocities,
    interface_crossing,
    normalized_temperature_rmse,
)


class TestStefanValidation(unittest.TestCase):
    def test_interface_crossing_interpolates_liquid_fraction(self):
        profile = [(0.0, 0.0), (0.4, 0.25), (0.8, 0.75), (1.0, 1.0)]
        self.assertAlmostEqual(interface_crossing(profile), 0.6)

    def test_exact_temperature_profile_has_zero_rmse(self):
        problem = CASES["matched"]
        time = 0.1
        state = interface_state(problem, time)
        profile = [
            (x, temperature(problem, x, time, state["lambda"]))
            for x in (0.0, state["position"], 0.01, 0.02)
        ]
        self.assertLess(
            normalized_temperature_rmse(problem, time, profile, state["lambda"]),
            1.0e-15,
        )

    def test_bulk_liquid_velocity_excludes_diffuse_interface_buffer(self):
        velocity = [(float(x), float(x)) for x in range(10)]
        liquid_fraction = [(float(x), 1.0) for x in range(10)]
        self.assertEqual(
            bulk_liquid_velocities(
                velocity, liquid_fraction, interface=1.5, interface_buffer_cells=3
            ),
            [5.0, 6.0, 7.0, 8.0, 9.0],
        )


if __name__ == "__main__":
    unittest.main()
