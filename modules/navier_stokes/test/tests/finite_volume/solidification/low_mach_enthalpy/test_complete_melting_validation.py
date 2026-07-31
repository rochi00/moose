import unittest

from validate_complete_melting import validate


class TestCompleteMeltingValidation(unittest.TestCase):
    def setUp(self):
        self.rows = [
            {
                "time": 0.5,
                "pcm_mass": 4.01,
                "pcm_liquid_volume": 0.3,
                "material_fraction": 0.58,
                "minimum_material_fraction": 0.0,
                "maximum_material_fraction": 1.0,
                "right_bulk_material_fraction": 0.0,
                "divergence_source_value": 0.2,
            },
            {
                "time": 1.0,
                "pcm_mass": 4.0,
                "pcm_liquid_volume": 0.66,
                "material_fraction": 2.0 / 3.0,
                "minimum_material_fraction": 0.0,
                "maximum_material_fraction": 1.0,
                "right_bulk_material_fraction": 0.0,
                "divergence_source_value": 0.0,
            },
        ]

    def test_conservative_complete_melting_result_passes(self):
        metrics = validate(
            self.rows,
            final_mass_relative_tolerance=0.005,
            transient_mass_relative_tolerance=0.015,
            volume_tolerance=0.01,
        )
        self.assertAlmostEqual(metrics["volume_error"], 0.0)

    def test_wrong_density_ratio_volume_fails(self):
        self.rows[-1]["material_fraction"] = 0.62
        with self.assertRaisesRegex(ValueError, "density-ratio"):
            validate(
                self.rows,
                final_mass_relative_tolerance=0.005,
                transient_mass_relative_tolerance=0.015,
                volume_tolerance=0.01,
            )

    def test_transient_mass_creation_fails(self):
        self.rows[0]["pcm_mass"] = 4.2
        with self.assertRaisesRegex(ValueError, "during melting"):
            validate(
                self.rows,
                final_mass_relative_tolerance=0.005,
                transient_mass_relative_tolerance=0.015,
                volume_tolerance=0.01,
            )


if __name__ == "__main__":
    unittest.main()
