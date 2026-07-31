#!/usr/bin/env python3

import argparse
import csv
import math


INITIAL_MATERIAL_FRACTION = 0.5
TIME_STEP = 1.0


def read_final_row(path):
    with open(path, newline="", encoding="utf-8") as csv_file:
        rows = list(csv.DictReader(csv_file))
    if not rows:
        raise ValueError(f"{path} contains no data rows")
    return {name: float(value) for name, value in rows[-1].items()}


def phase_fraction_balance(row):
    return (
        row["material_fraction"]
        - INITIAL_MATERIAL_FRACTION
        + TIME_STEP
        * (
            row["left_alpha_flux"]
            + row["right_alpha_flux"]
            - row["integrated_alpha_divergence"]
        )
    )


def validate(row, balance_tolerance, gas_bulk_tolerance):
    required_values = (
        "divergence_source_value",
        "integrated_alpha_divergence",
        "left_alpha_flux",
        "right_alpha_flux",
        "material_fraction",
        "minimum_material_fraction",
        "maximum_material_fraction",
        "left_bulk_material_fraction",
        "right_bulk_material_fraction",
    )
    if not all(math.isfinite(row[name]) for name in required_values):
        raise ValueError("VOF divergence validation contains a non-finite result")
    if row["divergence_source_value"] <= 0.0:
        raise ValueError("The manufactured case did not produce positive low-Mach divergence")
    if row["integrated_alpha_divergence"] <= 0.0:
        raise ValueError("The alpha-weighted low-Mach divergence source is not positive")
    if row["minimum_material_fraction"] < -1.0e-12:
        raise ValueError("The conservative VOF solve produced alpha below zero")
    if row["maximum_material_fraction"] > 1.0 + 1.0e-12:
        raise ValueError("The conservative VOF solve produced alpha above one")
    if abs(row["left_bulk_material_fraction"] - 1.0) > 1.0e-12:
        raise ValueError("The bulk PCM state was not preserved at alpha = 1")
    if row["right_bulk_material_fraction"] > gas_bulk_tolerance:
        raise ValueError(
            "The transported interface contaminated the far-gas bulk: "
            f"alpha = {row['right_bulk_material_fraction']:.8g}"
        )

    residual = phase_fraction_balance(row)
    if abs(residual) > balance_tolerance:
        raise ValueError(
            "The discrete phase-fraction balance did not close: "
            f"residual = {residual:.8g}"
        )
    return residual


def main():
    parser = argparse.ArgumentParser(
        description="Validate conservative VOF transport with nonzero low-Mach divergence."
    )
    parser.add_argument("--csv", default="vof_divergence_transport.csv")
    parser.add_argument("--balance-tolerance", type=float, default=1.0e-6)
    parser.add_argument("--gas-bulk-tolerance", type=float, default=5.0e-3)
    arguments = parser.parse_args()

    row = read_final_row(arguments.csv)
    residual = validate(
        row,
        balance_tolerance=arguments.balance_tolerance,
        gas_bulk_tolerance=arguments.gas_bulk_tolerance,
    )
    print(
        "VOF divergence validation: "
        f"alpha_min={row['minimum_material_fraction']:.8g}, "
        f"alpha_max={row['maximum_material_fraction']:.8g}, "
        f"far_gas_alpha={row['right_bulk_material_fraction']:.8g}, "
        f"balance_residual={residual:.8g}"
    )


if __name__ == "__main__":
    main()
