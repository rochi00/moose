import unittest

from validate_vof_divergence import phase_fraction_balance, validate


class TestVOFDivergenceValidation(unittest.TestCase):
    def setUp(self):
        self.row = {
            "divergence_source_value": 0.2,
            "integrated_alpha_divergence": 0.0525,
            "left_alpha_flux": 0.0,
            "right_alpha_flux": 0.0001,
            "material_fraction": 0.5524,
            "minimum_material_fraction": 0.0,
            "maximum_material_fraction": 1.0,
            "left_bulk_material_fraction": 1.0,
            "right_bulk_material_fraction": 0.001,
        }

    def test_phase_fraction_balance_uses_source_and_boundary_flux(self):
        self.assertAlmostEqual(phase_fraction_balance(self.row), 0.0)

    def test_valid_bounded_conservative_result_passes(self):
        self.assertAlmostEqual(
            validate(self.row, balance_tolerance=1.0e-12, gas_bulk_tolerance=0.005),
            0.0,
        )

    def test_unbounded_result_fails(self):
        self.row["maximum_material_fraction"] = 1.01
        with self.assertRaisesRegex(ValueError, "above one"):
            validate(self.row, balance_tolerance=1.0e-12, gas_bulk_tolerance=0.005)


if __name__ == "__main__":
    unittest.main()
