#!/usr/bin/env python3

"""Compare low-Mach Stefan CSV output with the analytical reference."""

import argparse
import csv
import glob
import math

from stefan_reference import CASES, interface_state, temperature


def read_profile(file_base, name, value_name):
    matches = sorted(glob.glob(f"{file_base}_{name}_*.csv"))
    if not matches:
        raise FileNotFoundError(f"no {name} CSV found for file base '{file_base}'")
    with open(matches[-1], newline="", encoding="utf-8") as stream:
        return [(float(row["x"]), float(row[value_name])) for row in csv.DictReader(stream)]


def interface_crossing(profile, level=0.5):
    for (x0, value0), (x1, value1) in zip(profile, profile[1:]):
        if value0 <= level <= value1 and value1 != value0:
            return x0 + (level - value0) * (x1 - x0) / (value1 - value0)
    raise ValueError(f"profile does not cross liquid fraction {level}")


def normalized_temperature_rmse(problem, time, profile, lam):
    error_squared = [
        (numerical - temperature(problem, x, time, lam)) ** 2 for x, numerical in profile
    ]
    scale = problem.temperature_initial - problem.temperature_cold
    return math.sqrt(sum(error_squared) / len(error_squared)) / scale


def bulk_liquid_velocities(
    velocity_profile, liquid_fraction_profile, interface, interface_buffer_cells=5
):
    """Return fully liquid samples outside the diffuse interface/drag transition."""
    if interface_buffer_cells < 0:
        raise ValueError("interface_buffer_cells must be nonnegative")
    spacings = [
        x1 - x0
        for (x0, _), (x1, _) in zip(velocity_profile, velocity_profile[1:])
        if x1 > x0
    ]
    if not spacings:
        raise ValueError("velocity profile must contain at least two distinct x coordinates")
    cutoff = interface + interface_buffer_cells * min(spacings)
    return [
        velocity
        for (x, velocity), (_, fraction) in zip(
            velocity_profile, liquid_fraction_profile
        )
        if x > cutoff and fraction > 0.999
    ]


def compare(case, time, file_base, velocity_interface_buffer_cells=5):
    problem = CASES[case]
    reference = interface_state(problem, time)
    liquid_fraction = read_profile(
        file_base, "liquid_fraction_profile", "liquid_fraction"
    )
    temperature_profile = read_profile(file_base, "temperature_profile", "T")
    velocity_profile = read_profile(file_base, "velocity_profile", "vel_x")

    numerical_interface = interface_crossing(liquid_fraction)
    interface_relative_error = abs(
        numerical_interface - reference["position"]
    ) / reference["position"]
    temperature_error = normalized_temperature_rmse(
        problem, time, temperature_profile, reference["lambda"]
    )

    liquid_velocities = bulk_liquid_velocities(
        velocity_profile,
        liquid_fraction,
        numerical_interface,
        velocity_interface_buffer_cells,
    )
    if not liquid_velocities:
        raise ValueError("no bulk-liquid velocity samples were found")
    velocity_error = max(
        abs(velocity - reference["liquid_velocity"]) for velocity in liquid_velocities
    )

    return {
        "interface": numerical_interface,
        "interface_reference": reference["position"],
        "interface_relative_error": interface_relative_error,
        "temperature_normalized_rmse": temperature_error,
        "liquid_velocity_reference": reference["liquid_velocity"],
        "liquid_velocity_max_error": velocity_error,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--file-base", required=True)
    parser.add_argument("--interface-relative-tolerance", type=float, required=True)
    parser.add_argument("--temperature-rmse-tolerance", type=float, required=True)
    parser.add_argument("--velocity-absolute-tolerance", type=float, required=True)
    parser.add_argument("--velocity-interface-buffer-cells", type=int, default=5)
    args = parser.parse_args()

    result = compare(
        args.case,
        args.time,
        args.file_base,
        args.velocity_interface_buffer_cells,
    )
    print(
        "interface={interface:.12g} reference={interface_reference:.12g} "
        "relative_error={interface_relative_error:.6g}".format(**result)
    )
    print(
        "temperature_normalized_rmse={temperature_normalized_rmse:.6g}".format(
            **result
        )
    )
    print(
        "liquid_velocity_reference={liquid_velocity_reference:.12g} "
        "max_error={liquid_velocity_max_error:.6g}".format(**result)
    )

    failures = []
    if result["interface_relative_error"] > args.interface_relative_tolerance:
        failures.append("interface-position tolerance")
    if result["temperature_normalized_rmse"] > args.temperature_rmse_tolerance:
        failures.append("temperature-RMSE tolerance")
    if result["liquid_velocity_max_error"] > args.velocity_absolute_tolerance:
        failures.append("liquid-velocity tolerance")
    if failures:
        parser.error("failed " + ", ".join(failures))


if __name__ == "__main__":
    main()
