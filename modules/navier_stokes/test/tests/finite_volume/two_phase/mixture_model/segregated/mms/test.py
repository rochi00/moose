import glob
import math
import os
import subprocess
import unittest

import pandas
from mooseutils import fuzzyAbsoluteEqual


# ---------------------------------------------------------------------------------------------
# Combined space and time refinement
#
# A pure time refinement on a fixed mesh plateaus once the temporal error drops below the spatial
# error, so the measured slope stops being the temporal order. Refining the mesh and the time step
# together with dt proportional to h avoids that: the total error is O(h^2) + O(dt^2), so a scheme
# that is second order in both gives a clean slope of two, while a first order time integrator pulls
# the combined slope down to one. The pair of cases below relies on exactly that contrast.
# ---------------------------------------------------------------------------------------------


# The all-terms study starts one level in. Its coarse end is not in the asymptotic range at all:
# between 10 and 20 the temperature converges at 1.39 and the velocities at 1.65, which drags the
# mean down without saying anything about the discretisation. From 20 the ladder is clean, the
# velocities running mildly super-convergent at 2.31 on the finest pair against 2.07 for the phase
# fraction and the temperature.
_FINE_LEVELS = [(20, 0.02), (40, 0.01), (80, 0.005)]


def _executable():
    candidates = glob.glob(os.path.join("..", "..", "..", "..", "..", "..", "..", "*-opt"))
    if not candidates:
        raise RuntimeError("could not locate the navier_stokes executable")
    return candidates[0]


_ALL_TERMS = "2d-drift-flux-mms-all-terms.i"
_ALL_FIELDS = ("L2phi", "L2u", "L2v", "L2T")

def _run_levels(extra_args, input_file, prefix, levels):
    """Run the transient case over the refinement levels, returning the final-time errors."""
    exe = _executable()
    rows = []
    for nx, dt in levels:
        base = "{}_{}".format(prefix, nx)
        cmd = [
            exe,
            "-i",
            input_file,
            "Mesh/gmg/nx={}".format(nx),
            "Mesh/gmg/ny={}".format(nx),
            "Executioner/dt={}".format(dt),
            "Outputs/file_base={}".format(base),
        ] + extra_args
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        rows.append(pandas.read_csv(base + ".csv").iloc[-1])
    return rows


def _order(rows, label):
    """Mean of the observed order over every consecutive pair of levels.

    Taking the finest pair alone makes the verdict hang on a single ratio, which is noisy: a study
    can be converging perfectly well and still put one pair outside a tolerance band, either
    because the coarse end has not reached the asymptotic range or because a field is briefly
    super-convergent. Averaging over the ladder uses every level that was paid for and is not
    thrown by one ratio.
    """
    orders = [
        math.log(rows[i - 1][label] / rows[i][label], 2.0) for i in range(1, len(rows))
    ]
    return sum(orders) / len(orders)


def check_temporal_convergence(
    extra_args,
    expected_order,
    tol=0.3,
    input_file=_ALL_TERMS,
    prefix="at",
    labels=_ALL_FIELDS,
    levels=_FINE_LEVELS,
):
    rows = _run_levels(extra_args, input_file, prefix, levels)
    for label in labels:
        value = _order(rows, label)
        print("%s, %f" % (label, value))
        assert fuzzyAbsoluteEqual(value, expected_order, tol)


class TestAllTermsBDF2(unittest.TestCase):
    """Every term of the model at once, with each transient in the conservative form the kernel
    assembles by default: d(rho_m u)/dt in momentum, d(rho_d alpha)/dt in the dispersed phase mass
    balance and d(rho_m cp_m T)/dt in energy, all with multipliers that move in time. Second order
    in every field, the energy included, which the studies above do not carry at all."""

    def test(self):
        check_temporal_convergence(
            ["Executioner/scheme=bdf2"],
            2.0,
            input_file=_ALL_TERMS,
            prefix="at_bdf2",
            labels=_ALL_FIELDS,
            levels=_FINE_LEVELS,
        )


class TestAllTermsImplicitEuler(unittest.TestCase):
    """The first order control for the case above. Without it a second order result would not
    establish that the study is sensitive to the time discretisation, which is the whole basis on
    which the transients are being verified."""

    def test(self):
        check_temporal_convergence(
            ["Executioner/scheme=implicit-euler"],
            1.0,
            input_file=_ALL_TERMS,
            prefix="at_ie",
            labels=_ALL_FIELDS,
            levels=_FINE_LEVELS,
        )
