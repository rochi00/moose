# Kokkos linear FV benchmark and parity harness

This harness supports the benchmark plan in `../docs/kokkos_fv_benchmark_plan.md`. It runs the
existing Kokkos linear FV advection-diffusion objects against their CPU counterparts and does
three things:

- checks **parity**: the CPU and Kokkos matrices and right-hand sides agree entrywise;
- checks **determinism**: repeated runs give bitwise-identical matrices and right-hand sides;
- measures **per-assembly timings**, broken down by section.

## Layout

| Path | Contents |
|---|---|
| `inputs/advdiff_kokkos.i` | 3D advection-diffusion-source on a generated mesh, Kokkos objects |
| `inputs/advdiff_cpu.i` | The same problem with the CPU linear FV objects, configured to match the Kokkos discretization |
| `inputs/advdiff_transient_{kokkos,cpu}.i` | Transient version (implicit Euler), Kokkos and CPU objects |
| `bench.py` | Driver with `parity`, `determinism`, `sweep` and `transient` subcommands |
| `results/` | JSON (parity, determinism) and CSV (sweep, transient, layout) outputs, time-stamped, plus result write-ups |
| `miniapp/` | Standalone Kokkos layout benchmark `fv_layout_bench` (variants V0-V4) and its sweep script |
| `meshes/make_meshes.sh` | Regenerates the mini-app meshes with the test object `FVMeshDumper` (about 2.1 GB, not committed) |

The design documents and dated findings are in `../docs/`.

## Building MOOSE for the harness

The harness needs MOOSE built from the `rc/kokkos_moose_navier_stokes` branch. That branch adds
two framework changes:

- **Timing sections** in `framework/src/kokkos/systems/KokkosLinearSystem.K`: `KokkosFVSyncToDevice`,
  `KokkosFVBoundaryData`, `KokkosFVElementalMatrix`, `KokkosFVFluxMatrix`, `KokkosFVElementalRHS`,
  `KokkosFVFluxRHS` and `KokkosFVSyncToHost`.
- **Thread override** in `framework/src/kokkos/base/KokkosMooseInit.K`: if `KOKKOS_NUM_THREADS` is
  set, Kokkos uses it instead of the MOOSE thread count. Linear FV variables require a single MOOSE
  thread, so without this override the OpenMP backend always runs one thread per rank.

The CPU (OpenMP) backend build used on macOS:

```bash
conda activate moose-next
cd ../..   # MOOSE repository root
./configure --with-kokkos=cpu
cd test && make -j14
```

On a GPU node, configure with `--with-kokkos` against a CUDA/HIP-enabled PETSc instead.

## Usage

```bash
conda activate moose-next
cd kokkos_fv/bench

# CPU vs Kokkos matrix/RHS parity (hex/tet/prism, orthogonal and skewed, 1 and 3 ranks)
python3 bench.py parity

# Bitwise reproducibility of the Kokkos assembly over 10 fresh runs, 4 OpenMP threads
python3 bench.py determinism --n 24 --threads 4 --runs 10

# Timing sweep: per-assembly time by section
python3 bench.py sweep --impl kokkos cpu --n 32 64 128 --threads 1 4 8 --np 1 --repeats 3
```

**Common options:**

| Option or variable | Effect |
|---|---|
| `--exe` or `BENCH_EXE` | MOOSE executable; default `../../test/moose_test-opt` |
| `--mpiexec` or `BENCH_MPIEXEC` | MPI launcher |
| `--workroot DIR` | where the temporary run directories go |
| `-v` | print every command |

## Method notes

**Repeated assembly.** The inputs use a `Transient` executioner with no time-derivative kernel.
Every step therefore reassembles and resolves the same steady system. The solve is one unmonitored
Richardson/Jacobi sweep (`l_max_its = 1`, `-ksp_norm_type none`), so it costs almost nothing and
always reports convergence.

**Per-assembly time.** Each case runs twice, with `--steps-low` and `--steps-high` steps. Section
times from `PerfGraphReporter` are differenced and divided by the step difference. This removes
setup and first-assembly costs such as COO preallocation.

**Inclusive times.** PerfGraph stores exclusive time per node; `bench.py` adds up each subtree. The
section common to both implementations is `LinearSystem::computeLinearSystemInternal`, reported as
`assembly`.

**Threads.** `--threads` sets `KOKKOS_NUM_THREADS` and `OMP_NUM_THREADS`. MOOSE always runs with
`--n-threads=1`. The CPU path ignores Kokkos threads.

**Parity tolerances.** Entrywise relative tolerance is `1e-12`. Absolute differences below
`1e-14 x` the largest entry are ignored. Mismatches are classified into rows of cells with a
boundary face and interior rows.

**Determinism.** Hashes the dumped PETSc matrix and RHS (`-ksp_view_mat` / `-ksp_view_rhs`
binary). Each run is a fresh process.

## Known discretization differences

These were found by `parity`.

- **Interior faces match** to round-off on every mesh tested: hex, tet and prism, orthogonal and
  skewed, serial and 3 ranks.
- **Boundary rows match only on orthogonal hex meshes.** Two causes:
  1. **Dirichlet distance.** The CPU `LinearFVAdvectionDiffusionFunctorDirichletBC` uses the
     normal projection `|d_Cf . n|` (`LinearFVBoundaryCondition::computeCellToFaceDistance`).
     `KokkosLinearFVFunctorDirichletBC` uses `|d_Cf|` (`datum.faceDCFMag()`).
  2. **Neumann face value.** The CPU `LinearFVAdvectionDiffusionFunctorNeumannBC` adds a
     tangential-gradient term `grad(phi_C) . (d_Cf - (d_Cf . n) n)` to the boundary value used
     by advection. The Kokkos Neumann BC returns the cell value.
- **Neumann meaning differs.** The Kokkos Neumann functor is a normal gradient; the CPU functor
  is a diffusive flux divided by `k`. The inputs use zero flux so the two coincide.
- **`KokkosLinearFVSource` default.** `scaling_factor` defaults to `1`, which builds a
  `KokkosConstantFunction`. The kernel rejects it because it requires a `KokkosParsedFunction`, so
  the input sets `scaling_factor` explicitly.

## Layout mini-app

```bash
conda activate moose-next
meshes/make_meshes.sh                  # needs moose_test-opt; writes meshes/*.fvd
cmake -S miniapp -B miniapp/build -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=$CONDA_PREFIX/lib/cmake/Kokkos
cmake --build miniapp/build
miniapp/run_layout_bench.sh            # writes results/layout_<timestamp>.csv
```
