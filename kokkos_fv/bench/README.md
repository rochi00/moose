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

# CPU vs Kokkos matrix/RHS parity (hex/tet/prism, orthogonal and skewed, nonorthogonal
# correction off and on, 1 and 3 ranks); --interp selects advected interpolation methods
python3 bench.py parity
python3 bench.py parity --interp upwind average minmod vanleer muscl

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

**Parity system.** `parity` compares the first assembled system. The steady inputs start from a
smooth non-uniform initial condition, so that the gradients, which the nonorthogonal correction and
the limited interpolation methods read, are nonzero in every cell, while both codes assemble from
the same solution. Comparing a later system would not work for the Minmod and Van Leer methods:
their flux-limiter ratio `r_f` divides by the difference of the cell values across a face, which
amplifies the round-off differences of the two earlier solves by orders of magnitude where that
difference is small.

**Parity tolerances.** Entrywise relative tolerance is `1e-12`. Absolute differences below
`1e-14 x` the largest entry are ignored. Mismatches are classified into rows of cells with a
boundary face and interior rows.

**Determinism.** Hashes the dumped PETSc matrix and RHS (`-ksp_view_mat` / `-ksp_view_rhs`
binary). Each run is a fresh process.

## Known discretization differences

These were found by `parity`.

- **The CPU and Kokkos systems match** to round-off on every mesh tested: hex, tet and prism,
  orthogonal and skewed, serial and 3 ranks (since phase 1 step 4, 2026-10-15), with and without
  the nonorthogonal correction (since step 5), and with each advected interpolation method (since
  step 6). Until then boundary
  rows differed on non-orthogonal meshes: the Kokkos Dirichlet BC used `|d_Cf|` instead of the
  normal projection `|d_Cf . n|`, and the Kokkos Neumann BC took a normal gradient and had no
  tangential-gradient term.
- **Gradient method.** The default Green-Gauss method computes the gradients of variables whose
  linear FV objects are all Kokkos objects with Kokkos, as only that sees the Kokkos boundary
  conditions, which the Neumann tangential term depends on.
- **Properties.** Since phase 1 step 7 the Kokkos objects accept numbers and
  `KokkosConstantFunction` objects for their properties, so `scaling_factor` no longer has to be
  set explicitly; the inputs still set it, as the CPU inputs do.

## Layout mini-app

```bash
conda activate moose-next
meshes/make_meshes.sh                  # needs moose_test-opt; writes meshes/*.fvd
cmake -S miniapp -B miniapp/build -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=$CONDA_PREFIX/lib/cmake/Kokkos
cmake --build miniapp/build
miniapp/run_layout_bench.sh            # writes results/layout_<timestamp>.csv
```

## Kokkos thread count and event counts

- `bench.py --threads` sets `KOKKOS_NUM_THREADS`. Outside the harness, `--kokkos-threads=<n>` on
  the MOOSE command line does the same and takes precedence. Linear FV problems cannot use
  `--n-threads` above 1 (see `../docs/findings_2026-10-14.md`, section 6).
- `tools/kokkos_event_counter.cpp` is a Kokkos Tools library that counts allocations, fences and
  kernel launches per calling MOOSE function:

  ```bash
  c++ -O2 -std=c++17 -shared -fPIC tools/kokkos_event_counter.cpp -o libkec.so
  KOKKOS_TOOLS_LIBS=$PWD/libkec.so KEC_OUT=counts.txt ../../test/moose_test-opt -i inputs/advdiff_kokkos.i n=8 steps=2
  ```

  Run with two step counts and subtract to get per-assembly counts.
