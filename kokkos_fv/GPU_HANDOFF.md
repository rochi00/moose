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

## Experiments

The first GPU campaign (V0-V4, transient with Kokkos types) is done. Its main conclusions shaped
the changes in `docs/findings_2026-10-14.md` and the expectations below: device-resident PETSc
types matter most, then a fixed host overhead of about 2 ms per assembly, then the coefficient
evaluation; V2 (atomics, precomputed slots) was the fastest assembly layout. The second campaign answers the open
design questions Q4, Q8, Q9 and Q11-Q19 (design document, section 15.2).

1. **Correctness first:** run `run_tests` as above and `python3 bench.py parity`. Expect the CPU
   results: orthogonal hex matches to round-off; on non-orthogonal meshes interior rows match and
   boundary rows differ for known reasons (README).
2. **Build the helpers** against the GPU PETSc/Kokkos:
   - the mini-app (`miniapp/build/fv_layout_bench`);
   - `tools/petsc_kokkos_probe` (build line in its header);
   - `tools/libkec.so` from `tools/kokkos_event_counter.cpp` (build line in its header and the
     README).
3. **Run the campaign:** `./run_gpu_campaign.sh [section ...]` from `kokkos_fv/bench`, inside a GPU
   allocation. Set `BENCH_MPIEXEC` (for example `srun`), `MESHES` (the mesh directory) and
   `NRANKS` (rank counts matching the node's GPUs, plus one oversubscribed count).

   | Section | Questions | What it runs |
   |---|---|---|
   | `miniapp` | Q11-Q14, Q18 | launch-cost microbenchmark; assembly layouts including the fused V2F; gradient variants G0, G0C, G1, G1C, G1CS, G2, G3 over four orderings |
   | `overhead` | Q19, Q9 | steady sweep with device and host PETSc types; Kokkos event counts per assembly; transient per-step times |
   | `types` | Q19 | checks the default `aijkokkos`/`kokkos` selection message and matrix type |
   | `q4` | Q4, Q9 | `-log_view` on 2-4 ranks; GPU-host copy counts of the ghost update and `MatSetValuesCOO` events |
   | `q8` | Q8 | the PETSc/Kokkos execution space probe |
   | `coefficients` | Q15 | parsed vs constant diffusion coefficient |
   | `solvers` | Q16 | Poisson: CG with Jacobi, BoomerAMG, GAMG with reuse; advection-diffusion: GMRES, BiCGStab, pipelined BiCGStab |
   | `ranks` | Q17 | rank scaling with device types and BiCGStab |

4. **Profile** where the section timers are not enough: Nsight Systems or rocprof on one
   `overhead` run, to separate launches, the constant-memory functor copies (Q7/Q14), and device
   execution (`KokkosFVWait`).

**Recording:** for each question, give the answer, the numbers, and the CSV or log it comes from,
in a dated findings file and in the design document's section 15.2 table.

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
