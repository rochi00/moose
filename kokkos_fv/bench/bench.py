#!/usr/bin/env python3
"""
Benchmark and parity harness for the Kokkos linear FV path in MOOSE.

Subcommands
  parity       Assemble the same advection-diffusion system with the CPU and the Kokkos objects
               and compare the matrices and right-hand sides entrywise.
  determinism  Assemble the Kokkos system repeatedly in fresh processes and check that the matrix
               and right-hand side are bitwise identical.
  sweep        Time the assembly for a range of mesh sizes, MPI ranks and Kokkos threads, and
               write per-assembly timings (by section) to a CSV file.

Per-assembly times use a difference method: each case runs with a low and a high step count, and
the section times are differenced, which removes setup and the first-assembly costs.

Run from any directory; paths are resolved relative to this file. Requires numpy.
"""

import argparse
import csv
import datetime
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
INPUTS = {
    "kokkos": os.path.join(HERE, "inputs", "advdiff_kokkos.i"),
    "cpu": os.path.join(HERE, "inputs", "advdiff_cpu.i"),
    "transient_kokkos": os.path.join(HERE, "inputs", "advdiff_transient_kokkos.i"),
    "transient_cpu": os.path.join(HERE, "inputs", "advdiff_transient_cpu.i"),
}
# Transient variants: (input key, extra overrides)
TRANSIENT_VARIANTS = {
    "kokkos": ("transient_kokkos", {}),
    "cpu": ("transient_cpu", {}),
    # Kokkos flux kernels with the CPU time derivative (host/device mixing)
    "mixed": (
        "transient_kokkos",
        {
            "LinearFVKernels/time/type": "LinearFVTimeDerivative",
            "Problem/allow_mixed_kokkos_linear_fv": "true",
        },
    ),
}
DEFAULT_EXE = os.path.normpath(os.path.join(HERE, "..", "..", "test", "moose_test-opt"))

# PETSc option prefix of the linear solve of system "u_sys" (determined from -options_left)
KSP_PREFIX = os.environ.get("BENCH_KSP_PREFIX", "")

# PerfGraph sections reported in sweeps (inclusive times). Sections absent for an implementation
# are left empty. "assembly" is the section common to the CPU and Kokkos paths.
SECTIONS = {
    "assembly": "LinearSystem::computeLinearSystemInternal",
    "cpu_kernels": "LinearSystem::LinearFVKernels_FullSystem",
    "cpu_gradients": "LinearSystem::LinearVariableFV_Gradients",
    "kokkos_total": "LinearSystem::computeKokkosLinearSystem",
    "kokkos_sync_to_device": "LinearSystem::KokkosFVSyncToDevice",
    "kokkos_bc_data": "LinearSystem::KokkosFVBoundaryData",
    "kokkos_face_coefficients": "LinearSystem::KokkosFVFaceCoefficients",
    "kokkos_elem_matrix": "LinearSystem::KokkosFVElementalMatrix",
    "kokkos_flux_matrix": "LinearSystem::KokkosFVFluxMatrix",
    "kokkos_elem_rhs": "LinearSystem::KokkosFVElementalRHS",
    "kokkos_flux_rhs": "LinearSystem::KokkosFVFluxRHS",
    "kokkos_wait": "LinearSystem::KokkosFVWait",
    "kokkos_sync_to_host": "LinearSystem::KokkosFVSyncToHost",
    "kokkos_mat_set_values": "KokkosSystem::MatSetValues",
    "linear_solve": "LinearSystem::LinearSystem::solve",
}


# ------------------------------------------------------------------------------------------------
# Running MOOSE
# ------------------------------------------------------------------------------------------------
def run_moose(args, impl, overrides, np_ranks, workdir, threads=None, extra=()):
    """Run one MOOSE case; returns (returncode, stdout+stderr)."""
    cmd = []
    if np_ranks > 1:
        cmd += [args.mpiexec, "-n", str(np_ranks)]
    cmd += [args.exe, "-i", INPUTS[impl], "--n-threads=1"]
    cmd += [f"{k}={v}" for k, v in overrides.items()]
    cmd += [f"Outputs/file_base={os.path.join(workdir, 'run')}"]
    cmd += list(extra)
    env = dict(os.environ)
    if threads is not None:
        env["KOKKOS_NUM_THREADS"] = str(threads)
        env["OMP_NUM_THREADS"] = str(threads)
    if args.verbose:
        print("  $", " ".join(cmd), f"(threads={threads})", flush=True)
    proc = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def check(rc, out, what):
    if rc != 0:
        sys.stderr.write(out[-4000:])
        raise SystemExit(f"MOOSE run failed: {what}")


# ------------------------------------------------------------------------------------------------
# PETSc binary format (big-endian; integer width detected from the header)
# ------------------------------------------------------------------------------------------------
MAT_FILE_CLASSID = 1211216
VEC_FILE_CLASSID = 1211214


def _int_dtype(raw, classid):
    if np.frombuffer(raw[:4], dtype=">i4")[0] == classid:
        return np.dtype(">i4")
    if np.frombuffer(raw[:8], dtype=">i8")[0] == classid:
        return np.dtype(">i8")
    raise ValueError("Not a PETSc binary object with the expected class id")


def read_petsc_mat(path):
    """Return a dense-free representation: dict {(row, col): value} plus the shape."""
    raw = open(path, "rb").read()
    it = _int_dtype(raw, MAT_FILE_CLASSID)
    off = 0

    def take(dtype, count):
        nonlocal off
        a = np.frombuffer(raw, dtype=dtype, count=count, offset=off)
        off += a.nbytes
        return a

    _, m, n, nz = take(it, 4)
    row_lengths = take(it, m)
    cols = take(it, nz)
    vals = take(np.dtype(">f8"), nz)
    rows = np.repeat(np.arange(m), row_lengths)
    return (int(m), int(n)), rows, cols.astype(np.int64), vals.astype(np.float64)


def read_petsc_vec(path):
    raw = open(path, "rb").read()
    it = _int_dtype(raw, VEC_FILE_CLASSID)
    n = np.frombuffer(raw, dtype=it, count=2)[1]
    return np.frombuffer(raw, dtype=">f8", count=int(n), offset=2 * it.itemsize).astype(np.float64)


def read_petsc_vecs(path):
    """Read every vector appended to a PETSc binary file (one per KSP solve)."""
    raw = open(path, "rb").read()
    it = _int_dtype(raw, VEC_FILE_CLASSID)
    out, off = [], 0
    while off < len(raw):
        cid, n = np.frombuffer(raw, dtype=it, count=2, offset=off)
        if cid != VEC_FILE_CLASSID:
            raise ValueError(f"unexpected class id {cid} at offset {off} in {path}")
        off += 2 * it.itemsize
        out.append(np.frombuffer(raw, dtype=">f8", count=int(n), offset=off).astype(np.float64))
        off += int(n) * 8
    return out


def mat_to_dict(rows, cols, vals, drop_zeros=True):
    d = {}
    for r, c, v in zip(rows.tolist(), cols.tolist(), vals.tolist()):
        if drop_zeros and v == 0.0:
            continue
        d[(r, c)] = d.get((r, c), 0.0) + v
    return d


def dump_options(workdir):
    p = KSP_PREFIX
    return [
        f"-{p}ksp_view_mat",
        f"binary:{os.path.join(workdir, 'A.bin')}",
        f"-{p}ksp_view_rhs",
        f"binary:{os.path.join(workdir, 'b.bin')}",
    ]


# ------------------------------------------------------------------------------------------------
# parity
# ------------------------------------------------------------------------------------------------
def compare_systems(dir_ref, dir_test, rtol, atol):
    """Compare A and b from two runs. Returns a dict of metrics."""
    shape_r, rr, cr, vr = read_petsc_mat(os.path.join(dir_ref, "A.bin"))
    shape_t, rt, ct, vt = read_petsc_mat(os.path.join(dir_test, "A.bin"))
    if shape_r != shape_t:
        return {"ok": False, "reason": f"shape {shape_r} vs {shape_t}"}
    a = mat_to_dict(rr, cr, vr)
    b = mat_to_dict(rt, ct, vt)
    keys = set(a) | set(b)
    # A row is a boundary row if it has fewer couplings than the most-connected row of its kind;
    # with a single element type this identifies cells with at least one boundary face.
    row_nnz = np.bincount(np.array([k[0] for k in a], dtype=np.int64), minlength=shape_r[0])
    boundary_row = row_nnz < row_nnz.max()
    max_abs, max_rel, n_bad = 0.0, 0.0, 0
    bad_rows = set()
    scale = max(abs(v) for v in a.values()) if a else 1.0
    for k in keys:
        x, y = a.get(k, 0.0), b.get(k, 0.0)
        diff = abs(x - y)
        rel = diff / max(abs(x), abs(y), 1e-300)
        max_abs = max(max_abs, diff)
        max_rel = max(max_rel, rel if diff > atol * scale else 0.0)
        if diff > atol * scale and rel > rtol:
            n_bad += 1
            bad_rows.add(k[0])
    rhs_r = read_petsc_vec(os.path.join(dir_ref, "b.bin"))
    rhs_t = read_petsc_vec(os.path.join(dir_test, "b.bin"))
    rhs_scale = max(np.max(np.abs(rhs_r)), 1e-300)
    rhs_diff = np.abs(rhs_r - rhs_t)
    rhs_rel = rhs_diff / np.maximum(np.maximum(np.abs(rhs_r), np.abs(rhs_t)), 1e-300)
    rhs_bad_mask = (rhs_diff > atol * rhs_scale) & (rhs_rel > rtol)
    rhs_bad = int(np.sum(rhs_bad_mask))
    mat_bad_rows = np.array(sorted(bad_rows), dtype=np.int64)
    return {
        "mat_bad_rows_boundary": int(np.sum(boundary_row[mat_bad_rows])) if len(mat_bad_rows) else 0,
        "mat_bad_rows_interior": int(np.sum(~boundary_row[mat_bad_rows])) if len(mat_bad_rows) else 0,
        "rhs_bad_boundary": int(np.sum(rhs_bad_mask & boundary_row)),
        "rhs_bad_interior": int(np.sum(rhs_bad_mask & ~boundary_row)),
        "ok": n_bad == 0 and rhs_bad == 0,
        "rows": shape_r[0],
        "nnz_ref": len(a),
        "nnz_test": len(b),
        "mat_max_abs": max_abs,
        "mat_max_rel": max_rel,
        "mat_bad": n_bad,
        "rhs_max_abs": float(np.max(rhs_diff)) if rhs_diff.size else 0.0,
        "rhs_bad": rhs_bad,
    }


def cmd_parity(args):
    results = []
    for elem in args.elem:
        for skew in args.skew:
            for n in args.n:
                for np_ranks in args.np:
                    case = f"{elem} n={n} skew={skew} np={np_ranks}"
                    dirs = {}
                    for impl in ("cpu", "kokkos"):
                        d = tempfile.mkdtemp(prefix=f"parity_{impl}_", dir=args.workroot)
                        ov = {"n": n, "elem_type": elem, "skew": skew, "steps": 1}
                        rc, out = run_moose(args, impl, ov, np_ranks, d, extra=dump_options(d))
                        check(rc, out, f"{impl} {case}")
                        if not os.path.exists(os.path.join(d, "A.bin")):
                            sys.stderr.write(out[-3000:])
                            raise SystemExit(
                                "Matrix dump not produced; set BENCH_KSP_PREFIX to the solver's "
                                "PETSc options prefix (see -options_left output above)."
                            )
                        dirs[impl] = d
                    m = compare_systems(dirs["cpu"], dirs["kokkos"], args.rtol, args.atol)
                    m.update({"case": case})
                    results.append(m)
                    status = "PASS" if m["ok"] else "FAIL"
                    print(
                        f"[{status}] {case}: rows={m.get('rows')} nnz cpu/kokkos="
                        f"{m.get('nnz_ref')}/{m.get('nnz_test')} "
                        f"mat max|d|={m.get('mat_max_abs', 0):.3e} rel={m.get('mat_max_rel', 0):.3e} "
                        f"bad={m.get('mat_bad')} (rows bnd/int {m.get('mat_bad_rows_boundary')}/"
                        f"{m.get('mat_bad_rows_interior')}) rhs max|d|={m.get('rhs_max_abs', 0):.3e} "
                        f"bad={m.get('rhs_bad')} (bnd/int {m.get('rhs_bad_boundary')}/"
                        f"{m.get('rhs_bad_interior')})",
                        flush=True,
                    )
                    if not args.keep:
                        for d in dirs.values():
                            shutil.rmtree(d, ignore_errors=True)
    write_json(args, "parity", results)
    return 0 if all(r["ok"] for r in results) else 1


# ------------------------------------------------------------------------------------------------
# determinism
# ------------------------------------------------------------------------------------------------
def cmd_determinism(args):
    hashes = []
    for i in range(args.runs):
        d = tempfile.mkdtemp(prefix="det_", dir=args.workroot)
        ov = {"n": args.n, "elem_type": args.elem, "steps": 1}
        rc, out = run_moose(args, "kokkos", ov, args.np, d, threads=args.threads, extra=dump_options(d))
        check(rc, out, f"determinism run {i}")
        h = hashlib.sha256()
        for f in ("A.bin", "b.bin"):
            h.update(open(os.path.join(d, f), "rb").read())
        hashes.append(h.hexdigest())
        shutil.rmtree(d, ignore_errors=True)
    unique = sorted(set(hashes))
    print(f"{args.runs} runs, {len(unique)} distinct matrix/RHS hashes")
    write_json(args, "determinism", {"hashes": hashes, "distinct": len(unique)})
    return 0 if len(unique) == 1 else 1


# ------------------------------------------------------------------------------------------------
# sweep
# ------------------------------------------------------------------------------------------------
def perf_sections(json_path):
    """
    Return {section_name: (inclusive_time_s, num_calls)} from a PerfGraphReporter JSON file.
    PerfGraph stores exclusive ("self") time per node; inclusive time adds the subtree. Sections
    appearing at several places in the graph are summed.
    """
    data = json.load(open(json_path))
    graph = None
    for step in reversed(data.get("time_steps", [])):
        g = step.get("perf", {}).get("graph")
        if g:
            graph = g
            break
    if graph is None:
        raise RuntimeError(f"No PerfGraphReporter graph in {json_path}")
    sections = {}

    def walk(name, node):
        """Accumulate and return the inclusive time of this node."""
        inclusive = float(node.get("time", 0.0))
        for k, v in node.items():
            if isinstance(v, dict):
                inclusive += walk(k, v)
        if "num_calls" in node:
            t, c = sections.get(name, (0.0, 0))
            sections[name] = (t + inclusive, c + int(node["num_calls"]))
        return inclusive

    for k, v in graph.items():
        if isinstance(v, dict):
            walk(k, v)
    return sections


def find_json(workdir):
    for f in os.listdir(workdir):
        if f.endswith(".json"):
            return os.path.join(workdir, f)
    raise RuntimeError(f"No JSON output in {workdir}")


def final_csv_row(workdir):
    """Last row of the CSV postprocessor output (cell count, memory)."""
    for f in os.listdir(workdir):
        if f.endswith(".csv"):
            with open(os.path.join(workdir, f)) as fh:
                rows = list(csv.DictReader(fh))
            return rows[-1] if rows else {}
    return {}


def cmd_sweep(args):
    rows = []
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    for impl in args.impl:
        for elem in args.elem:
            for n in args.n:
                for np_ranks in args.np:
                    for threads in args.threads:
                        per_rep = []
                        info = {}
                        for rep in range(args.repeats):
                            secs = {}
                            for steps in (args.steps_low, args.steps_high):
                                d = tempfile.mkdtemp(prefix=f"sweep_{impl}_", dir=args.workroot)
                                ov = {"n": n, "elem_type": elem, "steps": steps}
                                ov.update(parse_sets(args.set))
                                rc, out = run_moose(
                                    args, impl, ov, np_ranks, d, threads=threads,
                                    extra=args.petsc_args.split())
                                check(rc, out, f"sweep {impl} n={n} np={np_ranks} t={threads}")
                                secs[steps] = perf_sections(find_json(d))
                                info = final_csv_row(d)
                                shutil.rmtree(d, ignore_errors=True)
                            lo, hi = secs[args.steps_low], secs[args.steps_high]
                            diff = {}
                            for label, name in SECTIONS.items():
                                t_hi, c_hi = hi.get(name, (0.0, 0))
                                t_lo, c_lo = lo.get(name, (0.0, 0))
                                if c_hi > c_lo:
                                    diff[label] = (t_hi - t_lo) / (args.steps_high - args.steps_low)
                            per_rep.append(diff)
                        names = [k for k in SECTIONS if any(k in r for r in per_rep)]
                        row = {
                            "impl": impl,
                            "elem_type": elem,
                            "n": n,
                            "cells": int(float(info.get("ncells", 0))),
                            "np": np_ranks,
                            "threads": threads,
                            "repeats": args.repeats,
                            "max_mem_rank_mb": info.get("perf/max_memory_this_rank", ""),
                            "label": args.label,
                        }
                        for name in names:
                            vals = [r[name] for r in per_rep if name in r]
                            row[f"{name}_median_s"] = statistics.median(vals)
                            if len(vals) > 2:
                                q = statistics.quantiles(vals, n=4)
                                row[f"{name}_iqr_s"] = q[2] - q[0]
                        rows.append(row)
                        asm = row.get("assembly_median_s")
                        kk = row.get("kokkos_total_median_s")
                        print(
                            f"{impl:6s} {elem:6s} n={n:4d} cells={row['cells']:9d} np={np_ranks} "
                            f"t={threads}: assembly {asm * 1e3:9.3f} ms"
                            + (f" (kokkos part {kk * 1e3:9.3f} ms)" if kk else "")
                            + f"  = {asm / row['cells'] * 1e9:8.2f} ns/cell",
                            flush=True,
                        )
    os.makedirs(args.results, exist_ok=True)
    out_csv = os.path.join(args.results, f"sweep_{stamp}.csv")
    meta = ["label", "impl", "elem_type", "n", "cells", "np", "threads", "repeats", "max_mem_rank_mb"]
    timed = [f"{k}_{s}" for k in SECTIONS for s in ("median_s", "iqr_s")]
    present = {k for r in rows for k in r}
    fields = meta + [k for k in timed if k in present]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out_csv}")
    return 0


# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# transient
# ------------------------------------------------------------------------------------------------
def traj_diff(a, b):
    """Per-step max relative difference between two solution trajectories."""
    out = []
    for x, y in zip(a, b):
        scale = max(np.max(np.abs(x)), 1e-300)
        out.append(float(np.max(np.abs(x - y)) / scale))
    return out


def run_transient(args, variant, n, threads, np_ranks, workroot, tag):
    key, extra_ov = TRANSIENT_VARIANTS[variant]
    d = tempfile.mkdtemp(prefix=f"tr_{variant}_{tag}_", dir=workroot)
    ov = {"n": n, "elem_type": args.elem, "skew": args.skew, "steps": args.steps, "dt": args.dt}
    ov.update(extra_ov)
    ov.update(parse_sets(args.set))
    xfile = os.path.join(d, "x.bin")
    rc, out = run_moose(
        args,
        key,
        ov,
        np_ranks,
        d,
        threads=threads,
        extra=["-ksp_view_solution", f"binary:{xfile}::append", "-ksp_converged_reason"]
        + args.petsc_args.split(),
    )
    check(rc, out, f"transient {variant} n={n} t={threads} np={np_ranks}")
    its = [
        int(l.split("iterations")[-1].strip())
        for l in out.splitlines()
        if "Linear solve converged" in l and "iterations" in l
    ]
    secs = perf_sections(find_json(d))
    traj = read_petsc_vecs(xfile)
    info = final_csv_row(d)
    shutil.rmtree(d, ignore_errors=True)
    return secs, traj, its, info


def cmd_transient(args):
    rows = []
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    for n in args.n:
        for np_ranks in args.np:
            trajs = {}
            for variant in args.variants:
                for threads in args.threads if variant != "cpu" else [1]:
                    for rep in range(args.repeats if variant != "cpu" else 1):
                        secs, traj, its, info = run_transient(
                            args, variant, n, threads, np_ranks, args.workroot, f"r{rep}"
                        )
                        trajs[(variant, threads, rep)] = traj
                        steps = args.steps

                        def per_step(label):
                            t, c = secs.get(SECTIONS[label], (0.0, 0))
                            return t / steps if c else None

                        asm, ls = per_step("assembly"), per_step("linear_solve")
                        row = {
                            "variant": variant,
                            "elem_type": args.elem,
                            "skew": args.skew,
                            "n": n,
                            "cells": int(float(info.get("ncells", 0))),
                            "np": np_ranks,
                            "threads": threads,
                            "rep": rep,
                            "steps": steps,
                            "dt": args.dt,
                            "step_ms": secs.get("FEProblem::solve", (0.0, 0))[0] / steps * 1e3,
                            "assembly_ms": asm * 1e3 if asm else None,
                            "ksp_ms": (ls - asm) * 1e3 if (ls and asm) else None,
                            "kokkos_ms": (per_step("kokkos_total") or 0) * 1e3,
                            "host_fv_loop_ms": (per_step("cpu_kernels") or 0) * 1e3,
                            "host_gradients_ms": (per_step("cpu_gradients") or 0) * 1e3,
                            "ksp_its_mean": statistics.mean(its) if its else None,
                            "u_avg_final": info.get("u_avg"),
                            "petsc_args": args.petsc_args,
                            "label": args.label,
                        }
                        rows.append(row)
                        print(
                            f"{variant:6s} n={n} np={np_ranks} t={threads} rep={rep}: step "
                            f"{row['step_ms']:9.2f} ms  assembly {row['assembly_ms'] or 0:9.2f} ms  "
                            f"(kokkos {row['kokkos_ms']:8.2f}, host loop {row['host_fv_loop_ms']:8.2f}) "
                            f"ksp {row['ksp_ms'] or 0:9.2f} ms  its {row['ksp_its_mean']}",
                            flush=True,
                        )
            # trajectory comparisons against the CPU reference
            ref = trajs.get(("cpu", 1, 0))
            for key, traj in trajs.items():
                if key[0] == "cpu" or ref is None:
                    continue
                dv = traj_diff(ref, traj)
                print(
                    f"  {key[0]} t={key[1]} rep={key[2]} vs cpu: max rel diff first/last/max "
                    f"{dv[0]:.2e}/{dv[-1]:.2e}/{max(dv):.2e} over {len(dv)} steps"
                )
                rows.append(
                    {
                        "variant": f"{key[0]}_vs_cpu",
                        "n": n,
                        "np": np_ranks,
                        "threads": key[1],
                        "rep": key[2],
                        "traj_diff_first": dv[0],
                        "traj_diff_last": dv[-1],
                        "traj_diff_max": max(dv),
                    }
                )
            # run-to-run reproducibility of the non-CPU variants
            for (variant, threads, rep), traj in trajs.items():
                if rep == 0 or variant == "cpu":
                    continue
                dv = traj_diff(trajs[(variant, threads, 0)], traj)
                first = next((i for i, x in enumerate(dv) if x > 0), None)
                status = "identical" if first is None else f"differs from step {first + 1}"
                print(
                    f"  determinism {variant} t={threads}: rep {rep} vs rep 0 bitwise {status}; "
                    f"diff step1/last/max {dv[0]:.2e}/{dv[-1]:.2e}/{max(dv):.2e}"
                )
                rows.append(
                    {
                        "variant": f"{variant}_determinism",
                        "n": n,
                        "np": np_ranks,
                        "threads": threads,
                        "rep": rep,
                        "first_diff_step": (first + 1) if first is not None else "",
                        "traj_diff_first": dv[0],
                        "traj_diff_last": dv[-1],
                        "traj_diff_max": max(dv),
                    }
                )
    os.makedirs(args.results, exist_ok=True)
    out_csv = os.path.join(args.results, f"transient_{stamp}.csv")
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out_csv}")
    return 0


def parse_sets(sets):
    """Turn ["Functions/k_fn/expression=1", ...] into input-file overrides."""
    out = {}
    for item in sets or []:
        key, sep, value = item.partition("=")
        if not sep:
            raise SystemExit(f"--set expects key=value, got '{item}'")
        out[key] = value
    return out


def write_json(args, name, payload):
    os.makedirs(args.results, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    path = os.path.join(args.results, f"{name}_{stamp}.json")
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"Wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", default=os.environ.get("BENCH_EXE", DEFAULT_EXE))
    p.add_argument("--mpiexec", default=os.environ.get("BENCH_MPIEXEC", "mpiexec"))
    p.add_argument("--results", default=os.path.join(HERE, "results"))
    p.add_argument("--workroot", default=None, help="Directory for temporary run directories")
    p.add_argument("-v", "--verbose", action="store_true")
    sub = p.add_subparsers(dest="cmd", required=True)

    pp = sub.add_parser("parity")
    pp.add_argument("--n", type=int, nargs="+", default=[6, 10])
    pp.add_argument("--elem", nargs="+", default=["HEX8", "TET4", "PRISM6"])
    pp.add_argument("--skew", type=float, nargs="+", default=[0.0, 0.3])
    pp.add_argument("--np", type=int, nargs="+", default=[1, 3])
    pp.add_argument("--rtol", type=float, default=1e-12)
    pp.add_argument("--atol", type=float, default=1e-14, help="relative to the largest entry")
    pp.add_argument("--keep", action="store_true", help="keep the run directories")

    pd = sub.add_parser("determinism")
    pd.add_argument("--n", type=int, default=12)
    pd.add_argument("--elem", default="HEX8")
    pd.add_argument("--np", type=int, default=1)
    pd.add_argument("--threads", type=int, default=None)
    pd.add_argument("--runs", type=int, default=10)

    ps = sub.add_parser("sweep")
    ps.add_argument("--impl", nargs="+", default=["kokkos", "cpu"], choices=["kokkos", "cpu"])
    ps.add_argument("--elem", nargs="+", default=["HEX8"])
    ps.add_argument("--n", type=int, nargs="+", default=[16, 32, 64])
    ps.add_argument("--np", type=int, nargs="+", default=[1])
    ps.add_argument("--threads", type=int, nargs="+", default=[1])
    ps.add_argument("--repeats", type=int, default=3)
    ps.add_argument("--steps-low", type=int, default=3)
    ps.add_argument("--steps-high", type=int, default=13)
    ps.add_argument("--set", nargs="*", default=[], help="extra input overrides key=value")
    ps.add_argument("--petsc-args", default="", help="extra PETSc command-line options")
    ps.add_argument("--label", default="", help="free-form label stored in the CSV")

    pt = sub.add_parser("transient")
    pt.add_argument("--variants", nargs="+", default=["cpu", "kokkos", "mixed"],
                    choices=list(TRANSIENT_VARIANTS))
    pt.add_argument("--n", type=int, nargs="+", default=[32])
    pt.add_argument("--elem", default="HEX8")
    pt.add_argument("--skew", type=float, default=0.0)
    pt.add_argument("--np", type=int, nargs="+", default=[1])
    pt.add_argument("--threads", type=int, nargs="+", default=[1])
    pt.add_argument("--repeats", type=int, default=1, help="repeats of non-CPU variants (determinism)")
    pt.add_argument("--steps", type=int, default=50)
    pt.add_argument("--dt", type=float, default=0.01)
    pt.add_argument("--set", nargs="*", default=[], help="extra input overrides key=value")
    pt.add_argument("--petsc-args", default="", help="extra PETSc command-line options")
    pt.add_argument("--label", default="", help="free-form label stored in the CSV")

    args = p.parse_args()
    if not os.path.exists(args.exe):
        raise SystemExit(f"Executable not found: {args.exe} (use --exe or BENCH_EXE)")
    if args.workroot:
        os.makedirs(args.workroot, exist_ok=True)
    return {"parity": cmd_parity, "determinism": cmd_determinism, "sweep": cmd_sweep,
            "transient": cmd_transient}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
