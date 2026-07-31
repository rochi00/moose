#!/usr/bin/env python3

import argparse
import csv
import math


INITIAL_PCM_VOLUME = 0.5
SOLID_DENSITY = 8.0
LIQUID_DENSITY = 6.0
INITIAL_PCM_MASS = INITIAL_PCM_VOLUME * SOLID_DENSITY
EXPECTED_LIQUID_VOLUME = INITIAL_PCM_MASS / LIQUID_DENSITY


def read_rows(path):
    with open(path, newline="", encoding="utf-8") as csv_file:
        rows = list(csv.DictReader(csv_file))
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    return [{name: float(value) for name, value in row.items()} for row in rows]


def validate(
    rows,
    final_mass_relative_tolerance,
    transient_mass_relative_tolerance,
    volume_tolerance,
):
    completed_rows = [row for row in rows if row["time"] > 0.0]
    if not completed_rows:
        raise ValueError("The complete-melting result contains no completed time step")

    required_values = (
        "pcm_mass",
        "pcm_liquid_volume",
        "material_fraction",
        "minimum_material_fraction",
        "maximum_material_fraction",
        "right_bulk_material_fraction",
        "divergence_source_value",
    )
    if not all(
        math.isfinite(row[name]) for row in completed_rows for name in required_values
    ):
        raise ValueError("The complete-melting validation contains a non-finite result")

    if min(row["minimum_material_fraction"] for row in completed_rows) < -1.0e-12:
        raise ValueError("The coupled VOF solve produced alpha below zero")
    if max(row["maximum_material_fraction"] for row in completed_rows) > 1.0 + 1.0e-12:
        raise ValueError("The coupled VOF solve produced alpha above one")

    maximum_mass_error = max(
        abs(row["pcm_mass"] - INITIAL_PCM_MASS) for row in completed_rows
    )
    if maximum_mass_error > transient_mass_relative_tolerance * INITIAL_PCM_MASS:
        raise ValueError(
            "PCM mass conservation failed during melting: "
            f"maximum error = {maximum_mass_error:.8g}"
        )

    final = completed_rows[-1]
    final_mass_error = abs(final["pcm_mass"] - INITIAL_PCM_MASS)
    if final_mass_error > final_mass_relative_tolerance * INITIAL_PCM_MASS:
        raise ValueError(
            f"Final PCM mass conservation failed: error = {final_mass_error:.8g}"
        )
    if abs(final["material_fraction"] - EXPECTED_LIQUID_VOLUME) > volume_tolerance:
        raise ValueError(
            "The free surface did not reach the density-ratio volume prediction: "
            f"volume = {final['material_fraction']:.8g}"
        )

    liquid_volume_ratio = (
        final["pcm_liquid_volume"] / final["material_fraction"]
        if final["material_fraction"] > 0.0
        else 0.0
    )
    if liquid_volume_ratio < 0.97:
        raise ValueError(
            f"The PCM did not melt completely: liquid volume ratio = {liquid_volume_ratio:.8g}"
        )
    if abs(final["divergence_source_value"]) > 1.0e-12:
        raise ValueError("The final fully liquid state retained phase-change divergence")
    if final["right_bulk_material_fraction"] > 1.0e-3:
        raise ValueError("The moving interface contaminated the far-gas bulk")

    return {
        "final_mass_error": final_mass_error,
        "maximum_mass_error": maximum_mass_error,
        "volume_error": final["material_fraction"] - EXPECTED_LIQUID_VOLUME,
        "liquid_volume_ratio": liquid_volume_ratio,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Validate complete PCM melting and free-surface expansion."
    )
    parser.add_argument("--csv", default="complete_melting_free_surface.csv")
    parser.add_argument("--final-mass-relative-tolerance", type=float, default=5.0e-3)
    parser.add_argument(
        "--transient-mass-relative-tolerance", type=float, default=1.5e-2
    )
    parser.add_argument("--volume-tolerance", type=float, default=1.0e-2)
    arguments = parser.parse_args()

    metrics = validate(
        read_rows(arguments.csv),
        final_mass_relative_tolerance=arguments.final_mass_relative_tolerance,
        transient_mass_relative_tolerance=arguments.transient_mass_relative_tolerance,
        volume_tolerance=arguments.volume_tolerance,
    )
    print(
        "Complete-melting validation: "
        f"final_mass_error={metrics['final_mass_error']:.8g}, "
        f"maximum_mass_error={metrics['maximum_mass_error']:.8g}, "
        f"volume_error={metrics['volume_error']:.8g}, "
        f"liquid_volume_ratio={metrics['liquid_volume_ratio']:.8g}"
    )


if __name__ == "__main__":
    main()
