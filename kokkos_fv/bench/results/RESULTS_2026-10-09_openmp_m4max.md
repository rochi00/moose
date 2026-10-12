# Baseline results: Kokkos OpenMP backend, Apple M4 Max (2026-10-09)

**Setup.**
- Build: MOOSE `rc/kokkos_moose_navier_stokes` (from `next` at `410cc90c7e`, plus the timing sections and
  the `KOKKOS_NUM_THREADS` override), `--with-kokkos=cpu`, conda env `moose-next` (PETSc 3.25.4,
  Kokkos 4.7.04), opt build.
- Machine: Apple M4 Max (14 cores), 36 GB.
- Problem: 3D advection-diffusion-source, generated hex mesh, 1 MPI rank, Kokkos OpenMP threads as
  listed.
- Timings: per assembly, median of 3 repeats, difference method (steps 3 and 13).
- Raw data: `sweep_20261009-151625.csv` (Kokkos), `sweep_20261009-152830.csv` (CPU),
  `sweep_20261009-153347.csv` (constant-coefficient ablation).

These are CPU-backend numbers. They characterize the current code structure, not GPU performance.

## Parity (CPU objects vs Kokkos objects)

- **Orthogonal hex:** matrix and RHS identical to round-off (max |d| 4e-16), serial and 3 ranks.
- **Tet, prism and skewed hex:** all mismatches are in boundary rows. Interior rows match to
  round-off on every mesh. Two causes (details in `../README.md`):
  1. **Dirichlet distance:** `|d_Cf|` in Kokkos vs `|d_Cf . n|` on the CPU.
  2. **Neumann face value:** the CPU adds a tangential-gradient term; Kokkos uses the cell value.
- Raw data: `parity_20261009-151414.json`.

## Determinism (fresh process per run, 24^3 hex, 10 runs)

| OpenMP threads | Distinct matrix/RHS hashes |
|---|---|
| 1 | 1 |
| 4 | 10 |
| 8 | 10 |

The atomic accumulation in the current flux kernels makes threaded assembly non-reproducible.

## Per-assembly time (ms), 2,097,152 hex cells

| Section | Kokkos t=1 | Kokkos t=4 | Kokkos t=8 | CPU |
|---|---|---|---|---|
| **Total assembly** (`computeLinearSystemInternal`) | **1555** | **506** | **347** | **1874** |
| Kokkos flux matrix | 1252 | 325 | 176 | - |
| Kokkos elemental RHS (source) | 116 | 31 | 18 | - |
| Kokkos BC relation data | 14 | 5.4 | 4.4 | - |
| Kokkos sync to device / to host (incl. `MatSetValuesCOO`) | 3.0 / 12.4 | 2.1 / 10.0 | 2.0 / 9.6 | - |
| **Host linear FV loop with no CPU kernels attached** | **141** | **133** | **134** | - |
| CPU kernels (face + element loops) | - | - | - | 1552 |
| CPU gradients (`LinearVariableFV_Gradients`) | - | - | - | 628 |
| Peak memory per rank (MB) | 9937 | 10620 | 10256 | 5514 |

**Per cell:** Kokkos 741 / 241 / 166 ns per cell (t = 1 / 4 / 8); CPU 894 ns per cell. Per-cell
times are flat from 33k to 2.1M cells, as expected for a CPU backend.

## Constant-coefficient ablation (2,097,152 cells)

`k_fn = 1`, `s_fn = 1`; both are still `KokkosParsedFunction`s, but trivial ones.

| Section | t=1 parsed / constant | t=8 parsed / constant |
|---|---|---|
| Kokkos flux matrix | 1252 / 769 | 176 / 106 |
| Kokkos elemental RHS | 116 / 44 | 18 / 7 |
| Total assembly | 1555 / 984 | 347 / 268 |
| Host loop with no CPU kernels attached | 141 / 139 | 134 / 136 |

## Findings, mapped to the design

1. **The empty host linear FV loop is the largest fixed overhead.**
   - `ComputeLinearFVElementalThread` / `ComputeLinearFVFaceThread` iterate over every element and
     face even when every kernel is a Kokkos object. That costs about 135 ms per assembly at
     2.1M cells.
   - At 8 threads with constant coefficients, it costs **more than all the Kokkos work
     combined** (136 vs 132 ms).
   - It runs serially, so it limits thread scaling (Amdahl).
   - **Action:** skip the host loops when the system has no non-Kokkos linear FV objects. This is
     a small change in `LinearSystem::computeLinearSystemInternal`, and it supports R9 in the
     design.
2. **The flux-matrix pass dominates the Kokkos work** (80% at t=1), at about 60-100 ns per
   directed side.
   - Evaluating the parsed coefficient (device RPN interpreter) accounts for about 40% of that
     pass. The coefficient is evaluated 4 times per face: twice per directed side (diagonal and
     neighbor entry) and two directed sides per face.
   - **Supports:** D10 (precomputed cell/face property arrays), the compute-once face pass (D7),
     and the combined matrix + RHS hook (5.2 of the detailed design).
   - The remaining cost per directed side covers `Matrix::find` (binary search) twice, two
     atomics, and building the datum. It is the target of precomputed slots (H3) and the gather
     (H2), still to be measured.
3. **Thread scaling of the Kokkos passes is good.** The flux matrix speeds up 7.1x from 1 to 8
   threads. Total assembly speeds up only 4.5x because of the serial host loop (finding 1).
4. **Kokkos uses about 1.8x the memory of the CPU path** (about 2.1 KB per cell more at 2.1M
   cells). On the CPU backend the "device" copies are separate host allocations, so this number
   counts both host and device data. It should be broken down before GPU sizing (P12).
5. **The CPU path does more work than the Kokkos path.**
   - It spends 628 ms per assembly computing cell gradients. The CPU Neumann BC requests them
     (`LinearFVAdvectionDiffusionFunctorNeumannBC.C:34`) for its tangential correction, which the
     Kokkos Neumann BC does not have. That correction is parity difference 2.
   - So CPU vs Kokkos timings are not like-for-like. Treat the CPU column as a reference point,
     not a speedup baseline. Use the Kokkos-to-Kokkos comparisons (threads, ablations, future
     design switches) for decisions.
6. **The solve is negligible**, as intended: `linear_solve` minus assembly is under 4% everywhere.
