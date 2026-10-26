# GPU handoff: Kokkos linear FV benchmarks

Brief for a Claude Code session on a GPU machine. It checks the Kokkos linear FV implementation
(phases 0 and 1) on a GPU backend, then measures it and answers the open design questions. All
earlier results are from the Kokkos OpenMP backend on an Apple M4 Max, except the A100 mini-app
runs recorded in the design document.

## Repository and branch

- Fork `https://github.com/rochi00/moose.git`, branch `rc/kokkos_moose_navier_stokes` (based on
  `idaholab/moose` `next` at `410cc90c7e`). Add `idaholab/moose` as `origin` and the fork as `fork`.
- `CLAUDE.md` (a link to `AGENTS.md`) at the repo root holds the working rules.
- What the branch contains:
  - **Phase 0** (`docs/phase0_implementation_plan.md`): face-based assembly over
    `Moose::Kokkos::FVMesh` with face ranges, `FVSystemAssembly` (owned-row COO upload with
    `MatSetValuesCOO`, right-hand sides added into the PETSc vector arrays), atomic (default) and
    gather assembly (`Problem/kokkos_linear_fv_assembly`), one fence per assembly, and
    `aijkokkos`/`kokkos` PETSc types selected by default on GPU builds.
  - **Phase 1** (`docs/phase1_implementation_plan.md`): variables of other systems, the time
    derivative, `KokkosFVGreenGaussGradient` (with the Venkatakrishnan limiter), host-matching
    Dirichlet and Neumann BCs, non-orthogonal correction, the five advection interpolation methods,
    properties from numbers, functions, variables and `KokkosParsedFunctorMaterial`, and the
    Robin, outflow, scalar symmetry and extrapolated BCs.
- None of this has been built for CUDA or HIP yet.

## Read first

| File | For |
|---|---|
| `kokkos_fv/docs/kokkos_ns_linearfv_detailed_design.md` | The design; section 15.2 lists the open questions this campaign answers |
| `kokkos_fv/docs/phase0_implementation_plan.md`, `phase1_implementation_plan.md` | What was implemented and how it was verified on the CPU |
| `kokkos_fv/docs/findings_2026-10-14.md` | Latest CPU findings: per-assembly overhead, event counts, threading |
| `kokkos_fv/docs/kokkos_fv_benchmark_plan.md` | What to measure and why |
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
2. `cd test && make -j<N>`. This is the first device compile of the phase 0 and 1 code. Likely
   trouble spots: `RPNEvaluator::eval()` (its new trailing `arguments` pointer is shared with the
   Kokkos FE objects), the nested `acquireGradients()` with host staging of gradient vectors, and
   functor sizes (all linear FV objects take the constant-memory launch path, Q7).
3. Mini-app: `cmake -S kokkos_fv/bench/miniapp -B kokkos_fv/bench/miniapp/build
   -DKokkos_ROOT=<PETSc Kokkos install> -DCMAKE_BUILD_TYPE=Release`. With CUDA the compiler is
   normally Kokkos's `nvcc_wrapper`.
4. Helpers: `tools/petsc_kokkos_probe` and `tools/libkec.so` (build lines in the headers of
   `tools/petsc_kokkos_probe.cpp` and `tools/kokkos_event_counter.cpp`).
5. Meshes: `kokkos_fv/bench/meshes/make_meshes.sh <moose_test-opt>` (about 2.1 GB).

## Correctness

Do this before any timing; a wrong answer makes the timings meaningless.

1. **Regression tests**, from `test/`, serial and in parallel, then in gather mode:

   ```
   ./run_tests --re='kokkos' -j<N>
   ./run_tests --re='kokkos/(linearfv|fvmesh|fvgradientmethods)' -p 3 -j<N>
   ./run_tests --re='kokkos/(linearfv|fvmesh|fvgradientmethods)' -j<N> \
     --cli-args 'Problem/kokkos_linear_fv_assembly=gather'
   ```

   `-j` must be at least the `-p` count. The gold files are host results, so a pass means the
   device matches the CPU. Only `fvgradientmethods/green-gauss` uses the Kokkos gradient method
   (see the pitfalls); the other tests with gradients compute them on the host.
2. **Parity** (`run_gpu_campaign.sh correctness`): 240 cases comparing the Kokkos and host
   assembled systems entry by entry. On the CPU all 240 match, boundary rows included, so any
   mismatch on the GPU is a bug. The campaign stops if parity fails.
3. **Determinism** (`run_gpu_campaign.sh determinism`): atomic assembly is expected to give
   different systems from run to run; gather assembly must give one.

## Measurements

Run `./run_gpu_campaign.sh [section ...]` from `kokkos_fv/bench`, inside a GPU allocation. Set
`BENCH_MPIEXEC` (for example `srun`), `MESHES` (the mesh directory) and `NRANKS` (rank counts
matching the node's GPUs, plus one oversubscribed count). Without arguments it runs every section
below except `miniapp_gradients`.

| Section | Questions | What it runs |
|---|---|---|
| `correctness` | | parity, as above |
| `determinism` | | atomic and gather assembly, 5 runs each, skewed tet |
| `miniapp` | Q13, Q14 | launch-cost microbenchmark; assembly layouts including the fused V2F |
| `miniapp_gradients` | Q11, Q12, Q18 | gradient layouts over four orderings; answered on an A100, rerun only on a different GPU |
| `overhead` | Q19, Q9 | steady sweep with device and host PETSc types; Kokkos event counts per assembly; transient per-step times; steady sweep in gather mode |
| `phase1` | | skewed hex: non-orthogonal correction with host and Kokkos gradients, Van Leer (host gradient), MUSCL (limited Kokkos gradient), a diffusion coefficient from a `KokkosParsedFunctorMaterial` of `u` |
| `types` | Q19 | checks the default `aijkokkos`/`kokkos` selection message and matrix type |
| `q4` | Q4, Q9 | `-log_view` on 2-4 ranks; GPU-host copy counts of the ghost update and `MatSetValuesCOO` events |
| `q8` | Q8 | the PETSc/Kokkos execution space probe |
| `coefficients` | Q15 | parsed vs constant diffusion coefficient |
| `solvers` | Q16 | Poisson: CG with Jacobi, BoomerAMG, GAMG with reuse; advection-diffusion: GMRES, BiCGStab, pipelined BiCGStab |
| `ranks` | Q17 | rank scaling with device types and BiCGStab |

**Profile** one `overhead` run and one `phase1` run with Nsight Systems
(`nsys profile -t cuda,nvtx,mpi`, with the Kokkos Tools NVTX connector) or rocprof. Look for:

- mesh-sized host-device copies per assembly; there should be none with device PETSc types, except
  gradient staging when a gradient comes from a host method;
- the constant-memory functor copies before each launch (Q7, Q14);
- the gaps between launches, and how much of the assembly is `KokkosFVWait`.

**Recording:** for each question, give the answer, the numbers, and the CSV or log it comes from,
in a dated findings file and in the design document's section 15.2 table. For `phase1`, compare
the `cpu_gradients` and `kokkos_flux` columns across the labels.

## Known pitfalls

- **Timer semantics changed.** Linear FV kernels and BCs are dispatched without fences, and one
  fence (`KokkosFVWait`) precedes the upload. On a GPU, `KokkosFVBoundaryData`,
  `KokkosFVElemental` and `KokkosFVFlux` therefore measure launch time only, and the device work
  shows up in `KokkosFVWait`. The gradient computation fences itself, so `LinearVariableFV_Gradients`
  (bench column `cpu_gradients`, despite the name) is complete for host and Kokkos methods.
- **Gradients default to the host.** `MooseLinearVariableFV`'s `gradient_method` defaults to the
  host `green-gauss`. Select `KokkosFVGreenGaussGradient` in `[FVGradientMethods]` and set the
  variable's `gradient_method`, or on the command line:
  `FVGradientMethods/kgg/type=KokkosFVGreenGaussGradient Variables/u/gradient_method=kgg`.
  `FVAdvectedMinmodWeightBased` and `FVAdvectedVanLeerWeightBased` always request the host
  `green-gauss`, so their gradients are computed on the host and copied to the device every
  assembly. `FVAdvectedMUSCLDeferredCorrection` takes a `gradient_method`.
- MOOSE rejects an unprefixed `-mat_type`; use `-u_sys_mat_type aijkokkos` (system `u_sys`).
- Selecting only one of the device matrix and vector types crashes (a Kokkos matrix with host
  vectors). The default selection adds both, and nothing if the input or command line sets the
  system's matrix type or a vector type other than `kokkos`.
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

## Expectations

From the first GPU campaign and the A100 mini-app runs:

- Device-resident PETSc types matter most, then a fixed host overhead of about 2 ms per assembly
  (measured before the fence and allocation changes), then coefficient evaluation.
- The atomic face pass (V2, fused V2F) was 1.4-2.1x faster than the old per-side scheme; the
  gathers were slower. Gather mode should cost more on the GPU than on the CPU (1.2-1.3x there).
- Cell numbering changes assembly by at most 3% on good meshes; native numbering is kept.
- BiCGStab with Jacobi was 1.7x faster than GMRES(30) for advection-diffusion.

Not yet known: whether the remaining fences in `System::sync` and the constant-memory launch path
dominate at small sizes; the cost of the parsed-function evaluator on the device; and how much the
host gradients of Minmod and Van Leer cost per assembly.

## Recording results

Write results to `kokkos_fv/bench/results/` and a dated findings file in `kokkos_fv/docs/`, in the
same style as the CPU findings. Commit or push only when the user asks.
