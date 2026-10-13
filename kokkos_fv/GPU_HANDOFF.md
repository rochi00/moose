# GPU handoff: Kokkos linear FV benchmarks

Brief for a Claude Code session on a GPU machine. It repeats the CPU (Kokkos OpenMP) workflow from
`docs/findings_2026-10-09.md` and `docs/findings_2026-10-12.md` on a GPU backend.

## Repository and branch

- Fork `https://github.com/rochi00/moose.git`, branch `rc/kokkos_moose_navier_stokes` (based on
  `idaholab/moose` `next` at `410cc90c7e`). Add `idaholab/moose` as `origin` and the fork as `fork`.
- `CLAUDE.md` (a link to `AGENTS.md`) at the repo root holds the working rules.
- Commits on the branch:
  1. Speed up Kokkos linear FV assembly: skip the host linear FV loop when every linear FV
     object is Kokkos; error on mixed Kokkos/non-Kokkos linear FV systems (opt-out
     `Problem/allow_mixed_kokkos_linear_fv`); precomputed matrix slots; per-section timers;
     Kokkos thread count defaults to libMesh's.
  2. Test-only objects `FVMeshDumper` and `KokkosLinearFVImplicitEulerTest`.
  3. Design docs (`kokkos_fv/docs/`) and the harness (`kokkos_fv/bench/`).

## Read first

| File | For |
|---|---|
| `kokkos_fv/docs/kokkos_ns_linearfv_detailed_design.md` | The design and its hypotheses (H*, I*, R*) |
| `kokkos_fv/docs/kokkos_fv_benchmark_plan.md` | What to measure and why |
| `kokkos_fv/docs/findings_2026-10-12.md` | Latest CPU results, the implemented changes |
| `kokkos_fv/bench/results/LAYOUT_RESULTS_2026-10-09_openmp_m4max.md` | Mini-app variants V0-V4 on the CPU |
| `kokkos_fv/bench/README.md` | Harness usage |

## Questions to settle with the user before building

- Conda environment or module stack to use, and whether a GPU-enabled PETSc (with Kokkos and
  Kokkos-Kernels for CUDA or HIP) already exists or must be built.
- GPU model, scheduler, account, partition, and job time limits. Do not run GPU work on login
  nodes. Confirm job sizes before submitting.
- Work directory with about 20 GB free.

## Build

1. PETSc with Kokkos for the GPU (for example `--download-kokkos --download-kokkos-kernels
   --with-cuda` or `--with-hip`), then libMesh, then MOOSE configured with `--with-kokkos`
   (not `--with-kokkos=cpu`).
2. `cd test && make -j<N>`.
3. Correctness first: `./run_tests --re='kokkos|linearfv|linear_fv' -j<N>` and the Kokkos
   linear FV tests with `-p 3`. On the CPU build: 267 passed, 3 skipped, 0 failed.
4. Mini-app: `cmake -S kokkos_fv/bench/miniapp -B build -DKokkos_ROOT=<PETSc Kokkos install>
   -DCMAKE_BUILD_TYPE=Release`. It uses `Kokkos::DefaultExecutionSpace` and device Views
   throughout, so it should build for CUDA/HIP unchanged; with CUDA the compiler is normally
   Kokkos's `nvcc_wrapper`.
5. Meshes: `kokkos_fv/bench/meshes/make_meshes.sh <moose_test-opt>` (about 2.1 GB).

## Experiments, in order

1. **Parity:** `python3 bench.py parity`. Expect the CPU results: orthogonal hex matches to
   round-off; on non-orthogonal meshes interior rows match and boundary rows differ for known
   reasons (README).
2. **Steady sweep with section timers:** `python3 bench.py sweep --impl kokkos cpu --n 32 64 128
   --repeats 3`. Key question: how much of the Kokkos assembly is launch and fence overhead
   rather than compute. Compare section times against kernel times from a profiler (Nsight
   Systems, rocprof, or the Kokkos Tools kernel timer).
3. **Mini-app V0-V4:** `miniapp/run_layout_bench.sh`. Settles atomics (V2) against the gather
   (V3/V4), and per-face (V3) against half-face (V4) scratch, plus cell and face ordering.
4. **Transient with PETSc Kokkos types:** `python3 bench.py transient --variants cpu kokkos --n 64
   --repeats 3 --petsc-args "-vec_type kokkos -u_sys_mat_type aijkokkos -u_sys_ksp_type bcgs"`.
   Confirm the solve stays on the device (`-log_view` GPU flop share and copy counts).

## Known pitfalls

- MOOSE rejects an unprefixed `-mat_type`; use `-u_sys_mat_type aijkokkos` (system `u_sys`).
- `-ksp_view_solution` overwrites unless the viewer ends in `::append`.
- `KokkosLinearFVSource`'s default `scaling_factor` builds a `KokkosConstantFunction`, which the
  kernel rejects; the inputs set it explicitly.
- `preonly` KSP with a nonzero initial guess aborts; the inputs use Richardson + Jacobi with one
  iteration where only assembly is timed.
- PerfGraph JSON has exclusive times; `bench.py` sums them into inclusive times.
- Kokkos `Array::moveToDevice()` frees the host copy. The Kokkos mesh neighbor table is
  device-only, so host setup code must take neighbors from libMesh.
- Single runs are noisy; report medians of at least 3 repeats.
- `bench.py --threads` sets `KOKKOS_NUM_THREADS`/`OMP_NUM_THREADS`, which do not apply to a GPU
  backend; scale by MPI ranks (one GPU per rank) instead. Check how ranks are mapped to GPUs on
  the machine. `--mpiexec` can point at `srun` if the launcher accepts `-n`.
- MOOSE fences after each Kokkos dispatch, so PerfGraph section times include device work.

## Expectations from the CPU results (not yet measured on a GPU)

- Precomputed slots, skipping the host loop, and compute-once faces should help as much or more.
- Ordering still matters, through coalescing rather than cache reuse.
- Atomics vs gather and V3 vs V4 are open questions.
- Atomic assembly will likely be nondeterministic on every mesh; the gather is the deterministic
  option.
- Launch and fence overhead, the PETSc COO upload path, and parsed-function cost may dominate on
  a GPU and were invisible on the CPU.

## Recording results

Write results to `kokkos_fv/bench/results/` and a dated findings file in `kokkos_fv/docs/`, in the
same style as the CPU findings. Commit or push only when the user asks.
