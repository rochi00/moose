#!/usr/bin/env python3

"""Analytical reference for the variable-density two-phase Stefan problem."""

import argparse
import csv
import math
from dataclasses import dataclass


@dataclass(frozen=True)
class StefanProblem:
    rho_solid: float
    rho_liquid: float
    cp_solid: float = 910.0
    cp_liquid: float = 1042.4
    k_solid: float = 211.0
    k_liquid: float = 91.0
    latent_heat: float = 383840.0
    temperature_cold: float = 298.6
    temperature_initial: float = 973.6
    temperature_melt: float = 933.6
    temperature_reference: float = 933.6

    @property
    def density_ratio(self):
        return self.rho_solid / self.rho_liquid

    @property
    def diffusivity_solid(self):
        return self.k_solid / (self.rho_solid * self.cp_solid)

    @property
    def diffusivity_liquid(self):
        return self.k_liquid / (self.rho_liquid * self.cp_liquid)

    @property
    def effective_latent_heat(self):
        return self.latent_heat + (self.cp_liquid - self.cp_solid) * (
            self.temperature_melt - self.temperature_reference
        )


CASES = {
    "matched": StefanProblem(2475.0, 2475.0),
    "expansion": StefanProblem(500.0, 2700.0),
    "shrinkage": StefanProblem(2700.0, 500.0),
}


def stefan_residual(problem, lam, time, include_kinetic_energy=False):
    """Return left minus right sides of the interfacial energy balance."""
    if lam <= 0.0 or time <= 0.0:
        raise ValueError("lambda and time must be positive")

    alpha_s = problem.diffusivity_solid
    alpha_l = problem.diffusivity_liquid
    ratio = problem.density_ratio
    z_s = lam * math.sqrt(alpha_l / alpha_s)
    z_l = lam * ratio

    latent = problem.effective_latent_heat
    if include_kinetic_energy:
        latent -= 0.5 * (1.0 - ratio * ratio) * lam * lam * alpha_l / time

    left = problem.rho_solid * latent * lam * math.sqrt(alpha_l)
    solid_flux = (
        problem.k_solid
        * (problem.temperature_melt - problem.temperature_cold)
        * math.exp(-z_s * z_s)
        / (math.erf(z_s) * math.sqrt(math.pi * alpha_s))
    )
    liquid_flux = (
        problem.k_liquid
        * (problem.temperature_melt - problem.temperature_initial)
        * math.exp(-z_l * z_l)
        / (math.erfc(z_l) * math.sqrt(math.pi * alpha_l))
    )
    return left - solid_flux - liquid_flux


def solve_lambda(problem, time, include_kinetic_energy=False, tolerance=1.0e-13):
    """Solve for the first positive, physically continuous Stefan root."""
    low = 1.0e-10
    f_low = stefan_residual(problem, low, time, include_kinetic_energy)
    high = 0.25
    f_high = stefan_residual(problem, high, time, include_kinetic_energy)

    while f_low * f_high > 0.0 and high < 8.0:
        low, f_low = high, f_high
        high *= 2.0
        f_high = stefan_residual(problem, high, time, include_kinetic_energy)

    if f_low * f_high > 0.0:
        raise RuntimeError("failed to bracket the positive Stefan root")

    while high - low > tolerance * max(1.0, high):
        middle = 0.5 * (low + high)
        f_middle = stefan_residual(problem, middle, time, include_kinetic_energy)
        if f_low * f_middle <= 0.0:
            high = middle
        else:
            low, f_low = middle, f_middle

    return 0.5 * (low + high)


def interface_state(problem, time, include_kinetic_energy=False):
    lam = solve_lambda(problem, time, include_kinetic_energy)
    speed = lam * math.sqrt(problem.diffusivity_liquid / time)
    return {
        "lambda": lam,
        "position": 2.0 * lam * math.sqrt(problem.diffusivity_liquid * time),
        "speed": speed,
        "liquid_velocity": (1.0 - problem.density_ratio) * speed,
    }


def temperature(problem, x, time, lam=None, include_kinetic_energy=False):
    if x < 0.0:
        raise ValueError("x must be nonnegative")
    if lam is None:
        lam = solve_lambda(problem, time, include_kinetic_energy)

    interface = 2.0 * lam * math.sqrt(problem.diffusivity_liquid * time)
    if x <= interface:
        eta = x / (2.0 * math.sqrt(problem.diffusivity_solid * time))
        interface_eta = lam * math.sqrt(
            problem.diffusivity_liquid / problem.diffusivity_solid
        )
        return problem.temperature_cold + (
            problem.temperature_melt - problem.temperature_cold
        ) * math.erf(eta) / math.erf(interface_eta)

    eta = (
        x / (2.0 * math.sqrt(problem.diffusivity_liquid * time))
        - lam * (1.0 - problem.density_ratio)
    )
    return problem.temperature_initial + (
        problem.temperature_melt - problem.temperature_initial
    ) * math.erfc(eta) / math.erfc(lam * problem.density_ratio)


def write_profile(problem, time, length, points, output, include_kinetic_energy=False):
    state = interface_state(problem, time, include_kinetic_energy)
    with open(output, "w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("x", "temperature"))
        for index in range(points):
            x = length * index / (points - 1)
            writer.writerow(
                (f"{x:.16g}", f"{temperature(problem, x, time, state['lambda']):.16g}")
            )
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--length", type=float, default=1.0)
    parser.add_argument("--points", type=int, default=1281)
    parser.add_argument("--output", required=True)
    parser.add_argument("--include-kinetic-energy", action="store_true")
    args = parser.parse_args()

    if args.time <= 0.0:
        parser.error("--time must be positive")
    if args.length <= 0.0:
        parser.error("--length must be positive")
    if args.points < 2:
        parser.error("--points must be at least two")

    state = write_profile(
        CASES[args.case],
        args.time,
        args.length,
        args.points,
        args.output,
        args.include_kinetic_energy,
    )
    print(
        "lambda={lambda:.12g} interface={position:.12g} "
        "interface_speed={speed:.12g} liquid_velocity={liquid_velocity:.12g}".format(
            **state
        )
    )


if __name__ == "__main__":
    main()
