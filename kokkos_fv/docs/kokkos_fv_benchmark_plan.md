# Benchmark Plan: Validating the Kokkos Linear FV Design on Advection-Diffusion

**Goal.** Before porting any Navier-Stokes physics, measure whether the strategies in
`kokkos_ns_linearfv_detailed_design.md` actually pay off. The vehicle is the minimal Kokkos
linear FV implementation that already exists: `KokkosLinearFVAdvection`,
`KokkosLinearFVDiffusion`, `KokkosLinearFVSource`, and the functor Dirichlet/Neumann BCs.

This avoids two risks:

- building the whole NS stack on assumptions that turn out wrong;
- tuning code that later turns out not to matter.

**Baseline:** MOOSE `next` at `410cc90c7e`. Paths are relative to the `moose/` repository root.

---

## 1. What the minimal case can and cannot tell us

| Strategy (design reference) | Testable on advection-diffusion? | How |
|---|---|---|
| One fence per assembly instead of one per dispatch (R1, 4.3) | **Yes** | Directly |
| Face pass + gather vs the current directed-side atomic path (D7, R2, R3) | **Yes** | Directly for cheap kernels; with the proxy kernel (section 4) for NS-like heavy kernels |
| Precomputed matrix slots vs `Matrix::find` (4.1) | **Yes** | Directly |
| Owned-rows-only COO upload (4.1, D6) | **Yes** | Directly; needs MPI runs |
| Cell renumbering and face sorting (R4, 3.2) | **Yes** | Use meshes with scrambled numbering |
| Static BC memoization (5.3) | **Yes** | The Dirichlet and Neumann BCs here are static |
| Combined matrix + RHS hook (5.2) | **Yes** | Directly |
| Functor size and launch path (R10, Q7) | **Yes** | `sizeof` report plus launch overhead in Nsight |
| Device-resident PETSc types, no host round trips (5.3 of the overview) | **Yes** | Compare `-vec_type`/`-mat_type` host vs `kokkos` |
| Processor-face redundancy, owned-row assembly at scale (D5, D6) | **Yes** | Multi-GPU strong and weak scaling |
| AMG reuse/freeze policy for pressure (R8) | **Partly** | A diffusion-only system *is* a Poisson problem. Give it a coefficient that changes every step to imitate `rho * Ainv` |
| Gradient storage and Green-Gauss gather (6.2) | **Partly** | Only with the proxy gradient object (section 4) |
| Blocked halo exchanges (R6, Q4) | **Partly** | Standalone ghost-update microbenchmark (section 4.3) |
| Rhie-Chow, momentum sharing (R7), SIMPLE sync count | **No** | Wait for phase 2 |

---

## 2. Hypotheses and decision rules

Each strategy is stated as a prediction with a numeric threshold, so the result decides whether it
is adopted.

"Small" means at most 256k cells per GPU and "large" means at least 4M cells per GPU. Times are
medians over repeated assemblies (section 5.4).

| ID | Hypothesis | Expected outcome | Adopt if |
|---|---|---|---|
| H1 | Removing per-dispatch fences cuts assembly time where launches dominate | 20-50% less assembly time at small sizes; under 5% change at large sizes | At least 15% gain at small sizes, and no regression above 3% elsewhere |
| H2a | Face pass + gather is no slower than directed-side atomics for **cheap** kernels | Within ±15% at large sizes | Not more than 15% slower (we accept that cost for determinism) |
| H2b | Face pass + gather beats directed-side atomics for **heavy** face kernels (proxy) | Faster once face work exceeds about 2x the cheap kernels' | Crossover found at or below NS-like intensity (section 4.1) |
| H2c | Gather mode is bitwise deterministic; atomic modes are not | Matrix hash identical over 10 runs (gather); differs (atomics, on GPU) | Required for gather; just recorded for atomics |
| H3 | Precomputed slots beat `Matrix::find` | 10-30% less face-pass time | At least 10% |
| H4 | Owned-rows-only COO reduces upload time and memory | `MatSetValuesCOO` time and device memory down in multi-rank runs | Any reduction, with no correctness change |
| H5 | Renumbering helps scrambled meshes and does not hurt ordered ones | 1.1-1.7x on scrambled meshes; within ±3% on generated meshes | At least 1.1x on scrambled meshes and no regression above 3% |
| H6 | BC memoization removes BC launches | BC time goes to about 0 after the first assembly | Always (correctness check only) |
| H7 | Kokkos objects are small enough for the fast launch path | `sizeof` below the parameter-passing limit of the toolchain | If not: escalate (design 5.7) |
| H8 | Device PETSc types eliminate mesh-sized host traffic | Nsight shows 0 bytes of mesh-sized H2D/D2H per assembly | Required |
| H9 | Strong scaling holds to about 1-2M cells per GPU | At least 70% parallel efficiency of *assembly* down to about 1M cells per GPU | Recorded, and used to set the target sizes |
| H10 | AMG freeze/reuse amortizes setup on a Poisson problem with a changing coefficient | Total solve time down at least 30% at fixed accuracy | Policy parameters tuned from the results |

---

## 3. Implementation for the experiments

### 3.1 Switches, not branches

Every strategy is a runtime switch (or a compile-time option where unavoidable) in **one build**.
The baseline is then the same binary with all switches off. Proposed parameters, all
experimental and hidden from normal users:

| Parameter (where) | Values | Default for the experiments |
|---|---|---|
| `kokkos_fv_assembly` (system or Problem) | `directed_atomic` (current code) / `face_atomic` / `face_gather` | `directed_atomic` |
| `kokkos_dispatch_fence` (Problem) | `per_dispatch` (current) / `per_assembly` | `per_dispatch` |
| `kokkos_fv_matrix_slots` (system) | `search` (current) / `precomputed` | `search` |
| `kokkos_fv_coo_rows` (system) | `owned_and_ghost` (current) / `owned` | `owned_and_ghost` |
| `kokkos_cell_ordering` (Mesh) | `native` / `morton` / `rcm` | `native` |
| `kokkos_fv_bc_memoize` (BC) | `false` / `true` | `false` |
| `kokkos_fv_combined_hooks` (kernel) | `false` / `true` | `false` |

- **Parity first.** Each switch value must pass the correctness gate (section 5.1) before it is
  timed.
- **Leave current code intact.** `directed_atomic` must remain the unchanged current code path, so
  it is a true baseline.

### 3.2 Finer timing sections

Today the whole Kokkos linear assembly is one section,
`TIME_SECTION("computeKokkosLinearSystem", 1)` (`framework/src/kokkos/systems/KokkosLinearSystem.K:101`).
Add sub-sections:

| Section | Covers |
|---|---|
| `KokkosFVSyncToDevice` | `sync(HOST_TO_DEVICE)` and `copyToDevice` |
| `KokkosFVBC` | BC relation passes |
| `KokkosFVFace` | flux kernel passes |
| `KokkosFVCell` | elemental kernel passes |
| `KokkosFVGather` | gather (new path only) |
| `KokkosFVUpload` | `sync(DEVICE_TO_HOST)`, including `MatSetValuesCOO` and vector close |

**Accuracy caveat.** With fences removed, a section only measures *launch* time, not GPU time. So:

- for **per-pass** timing, run with `kokkos_dispatch_fence = per_dispatch`, or put a fence at each
  section end under a `--kokkos-timing-fences` flag;
- for **whole-assembly** comparisons, fence once at the end.

Report both.

### 3.3 Counters

Kokkos Tools `space-time-stack` gives kernels and fences per assembly. Add a counter printout for:

- kernel launches;
- fences;
- bytes copied host to device and device to host by the Kokkos `Array` copy APIs.

Pass criterion for H8: the byte counts stay O(1) per assembly, not O(cells).

### 3.4 Functor size report (H7)

Under a debug flag, `DispatcherRegistry::addDispatcher` prints `sizeof(Object)` for each
registered Kokkos type. Compare against the launch-path thresholds of the Kokkos build (design
5.7), and confirm with Nsight Systems: on the constant-memory path, look for memcpy-to-symbol and
event-sync calls before each launch.

---

## 4. Proxy objects for NS-like workload (test-only, under `test/src/kokkos/`)

The production advection/diffusion kernels are extremely light: a few reads and a parsed function
per face. NS momentum faces read many cell fields. These proxies imitate that without porting NS.

### 4.1 `KokkosLinearFVProxyFlux`

A flux kernel with two tunable parameters:

- `n_cell_fields`: how many extra cell arrays to read on each side of the face (emulates velocity
  components and gradients). Each is an `Array<Real>` of size `n_owned + n_ghost`.
- `flop_factor`: a loop of dependent FMAs per face (emulates interpolation and NOC arithmetic).

It produces a valid, diagonally dominant contribution (diffusion-like plus a small term derived
from the extra fields), so the system remains solvable and results are checkable.

**Calibration.** The NS momentum flux with deviatoric terms reads about:

- per side: 3 velocities + 9 gradient components + mu;
- per face: `mdot`.

That corresponds to `n_cell_fields` of about 13 with `flop_factor` of about 50-100 FMAs per face.
Sweep `n_cell_fields` over {0, 4, 13, 24} and `flop_factor` over {0, 50, 200} to find where the
face pass + gather crossover lies (H2b).

### 4.2 `KokkosLinearFVProxyGradient`

A cell gather over `elem_faces` computing a Green-Gauss gradient of the solution into a block
vector. This measures:

- the gather's memory behavior with and without renumbering (H5);
- the time of the gradient pass that every NS equation will need.

### 4.3 Ghost-update microbenchmark (Q4, R6)

A small test object or standalone PETSc program on the target machine. It times:

- `VecGhostUpdateBegin/End` for `k` separate ghosted vectors against one block vector with
  `bs = k`, for `k` in {1, 3, 9, 13};
- host vector types (`VECMPI`) against Kokkos vector types (`VECMPIKOKKOS` with `VecMPISetGhost`),
  checking in Nsight whether data stays on the device.

This answers Q4 and quantifies R6 without any assembly code.

---

## 4A. Indexing and layout experiments

The indexing strategies in the detailed design (sections 3-4 there) are decisions about data
layout. Each one has realistic alternatives, and they interact: for example, the best face order
depends on the cell numbering. They are tested separately from the assembly-algorithm switches
in section 3.1.

### 4A.1 Decisions under test

| ID | Decision (design ref.) | Proposed | Alternatives to compare | Expected effect |
|---|---|---|---|---|
| I1 | Face storage (D1, 3.4) | One entry per face, SoA | Current per-(side, cell) directed storage with every interior face twice | About half the geometry memory. Face pass faster for heavy kernels, same for light ones |
| I2 | Face→cell references (D2) | `elem[f]`, `neighbor[f]` as `CellID` | Current `neighbor(side, c)` lookup; 32-bit vs 64-bit index (`--with-kokkos-index-size=4|8`) | 32-bit saves bytes per face; recorded, no decision expected |
| I3 | Cell→face table (D3) | Padded `(side, c)`, LEFT layout (a cell's sides contiguous) | (a) padded with RIGHT layout (one side across all cells contiguous); (b) CSR `offset[c]` + `faces[k]`; (c) per-block stride | Padded LEFT is best for hex/tet; CSR is better only for mixed topologies |
| I4 | Orientation inside the gather (3.3) | Implicit `elem[f] == c` | (a) stored `int8` sign per half-face; (b) sign packed into the `FaceID` high bit | Small. Matters only if the gather is load-bound |
| I5 | Assembly scratch layout (R2, 4.4) | Half-face scratch `(side, c)`, so the gather reads contiguously | Per-face scratch `diag_E, diag_N, rhs_E, rhs_N` (original D7), where the gather reads scattered | Half-face scratch reduces gather DRAM traffic and transactions |
| I6 | Face order within ranges (3.5) | Sorted by owner cell | (a) `FaceInfo` creation order; (b) sorted by min(cell IDs); (c) Morton key of the face centroid | Owner sort is best together with cell renumbering |
| I7 | Cell numbering (R4, 3.2) | Morton (or RCM) within each block | native libMesh order | 1.1-1.7x on poorly ordered meshes |
| I8 | Worklist representation (3.6) | Range list + binary-search dispatch | (a) one launch per range; (b) explicit `FaceID` worklist array (today's approach, but face-indexed) | Range list wins when a variable has many ranges (many blocks or sidesets); otherwise a tie |
| I9 | Geometry representation (3.4) | `Real3` per face plus precomputed `inv_dn` and `noc` | (a) separate x/y/z arrays; (b) recomputing `dCN`, `inv_dn` and `noc` from centroids in the kernel | Precomputing saves arithmetic but adds bytes read per face. It wins when the face pass is compute-bound and loses when it is bandwidth-bound. Must be measured |
| I10 | COO slot order (4.1) | Row-major, gather writes each row contiguously | Face-order slots (the face pass writes off-diagonals directly) | Row-major gives contiguous writes and a trivial PETSc COO path |

### 4A.2 Two-level approach

Changing a core layout inside MOOSE for each variant is slow to build and couples every variant
to the framework. So:

1. **Standalone layout mini-app** (Kokkos only, no MOOSE or PETSc in the timed loop).
   - It reads the mesh index and geometry arrays from files dumped by MOOSE.
   - It implements the face pass, the gather, and a Green-Gauss gradient in each layout variant
     (I3, I4, I5, I6, I7, I9, I10, and I8 where it applies).
   - Compiles in seconds, so a full factorial is cheap.
   - Used to pick the winners.
2. **In-MOOSE confirmation** of only the winners and the decisions that need the framework:
   - I1 (directed vs face storage, which is the H2 comparison);
   - I7 (renumbering across FE and FV);
   - I8 (block-restricted variables and real BC ranges);
   - I10 (the actual `MatSetValuesCOO` path).

**Mesh dump.** A test-only MOOSE object writes per-rank binary files from `FaceInfo`/`ElemInfo`:
cell volumes and centroids, face cells, sides, normals, areas, `gC`, `dCN`, block and boundary
IDs, owned/ghost flags. Run it with `-p n` to get real partitions, including processor faces. The
mini-app and MOOSE then use identical meshes.

**Mesh set:**

- generated hex;
- tet (from a generated hex split, or an imported tet mesh);
- mixed hex/prism;
- skewed hex;
- a **scrambled** hex, the same geometry with element IDs randomly permuted, which is the worst
  case for I6/I7;
- a mesh with 1, 8 and 64 blocks and many sidesets, for I8.

### 4A.3 Measurements per variant

| Metric | Tool | Why |
|---|---|---|
| Time of face pass, gather and gradient pass (median of at least 50 repetitions, fenced) | Kokkos timer in the mini-app; MOOSE timing sections for confirmation | Primary |
| DRAM bytes read/written per pass (`dram__bytes_read.sum`, `dram__bytes_write.sum`) | Nsight Compute | Separates "fewer bytes" from "better access pattern" |
| Global load efficiency and sectors per request | Nsight Compute | Measures coalescing (I3, I5, I6, I9) |
| L2 hit rate | Nsight Compute | Locality from renumbering and face order (I6, I7) |
| Achieved bandwidth as a fraction of STREAM | derived | Comparable across variants and GPUs |
| Device memory for the index and geometry arrays | allocation sizes | I1, I2, I3, I9 |
| Index build (setup) time on the host | timer | Renumbering and sorting cost must stay small next to one solve |

**Byte model.** For each variant, write down the theoretical bytes per face and per cell. If
measured DRAM bytes are far above the model, the access pattern is wasting cache lines. That is a
diagnosis, not just a timing.

### 4A.4 Correctness invariants (checked for every variant)

- Every interior face appears exactly once in the face arrays and exactly twice in the cell→face
  table (once per owned side).
- The signed area vectors of each closed cell sum to about 0.
- `gC`, `dCN`, area and normal match `FaceInfo` exactly (they are copied, not recomputed), except
  in variant I9(b), where the tolerance is round-off.
- Face pass + gather produces the same matrix and RHS in every layout variant: bitwise for
  variants that keep the summation order, round-off otherwise.
- Renumbering is a bijection and preserves block ranges.

### 4A.5 Order and decision rule

1. Mini-app on one GPU, single rank, all meshes, at 1M and 16M cells:
   - full factorial of I6 x I7 (3 x 3);
   - then I3, I4, I5, I9 and I10 one at a time on the best I6/I7 combination.
2. Adopt a variant if it is at least 5% faster on at least one mesh, with no slowdown above 3% on
   the others. Alternatively, adopt it if it is equally fast and uses at least 20% less memory.
3. Confirm the chosen combination in MOOSE (I1, I7, I8, I10) with the H2/H5 runs of section 5.

## 5. Methodology

### 5.1 Correctness gate (run before every timing campaign)

1. **Matrix and RHS parity.** For each switch combination, assemble once and compare the PETSc
   matrix and RHS entrywise against:
   - the `directed_atomic` baseline;
   - the CPU objects (`LinearFVAdvection`, `LinearFVDiffusion`, `LinearFVSource`, and the CPU
     functor BCs, as in `test/tests/kokkos/linearfvkernels/diffusion/match-non-kokkos.i`).

   Use relative tolerance 1e-12 on each nonzero, ignoring explicit zeros.
2. **Solution parity.** The existing tests
   `test/tests/kokkos/linearfvkernels/{diffusion,advection-diffusion}` pass with each switch
   combination, with `-p 1` and `-p 3`.
3. **Determinism (H2c).** Assemble 10 times in a fresh process each time; hash the COO value array
   and the RHS (for example a 64-bit FNV over the raw bytes).

### 5.2 Benchmark inputs

A family of 3D steady advection-diffusion problems built from the existing objects:

```
[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${n}  ny = ${n}  nz = ${n}
    elem_type = HEX8               # also TET4 / PRISM6 variants
  []
  # variants: a perturbed (skewed) mesh via a mesh generator that jitters node positions,
  # and a "scrambled" mesh: same geometry, element IDs randomly permuted before
  # Kokkos numbering (e.g. via a test-only mesh generator)
[]

[Problem]
  linear_sys_names = 'u_sys'
[]

[Variables]
  [u]
    type = MooseLinearVariableFVReal
    solver_sys = 'u_sys'
  []
[]

[LinearFVKernels]
  [adv]
    type = KokkosLinearFVAdvection
    variable = u
    velocity = '1 0.5 0.25'
  []
  [diff]
    type = KokkosLinearFVDiffusion
    variable = u
    diffusion_coeff = k_fn
  []
  [src]
    type = KokkosLinearFVSource
    variable = u
    source_density = s_fn
  []
  # proxy variant: replace adv/diff with KokkosLinearFVProxyFlux
[]

[LinearFVBCs]
  [in]
    type = KokkosLinearFVFunctorDirichletBC
    variable = u
    boundary = 'left bottom back'
    functor = g_fn
  []
  [out]
    type = KokkosLinearFVFunctorNeumannBC
    variable = u
    boundary = 'right top front'
    functor = zero_fn
  []
[]

[Executioner]
  type = Transient                 # repeated assembly; each step reassembles
  num_steps = ${steps}             # e.g. 60; discard the first 10 as warm-up
  dt = 1
  system_names = u_sys
  # Assembly campaign: make the solve negligible
  petsc_options_iname = '-ksp_type -pc_type -ksp_max_it'
  petsc_options_value = 'richardson jacobi 1'
  # End-to-end campaign: real solver
  # petsc_options_value = 'bcgs jacobi 200'   (advection-diffusion)
  # Poisson campaign (H10): diffusion only, k_fn time dependent, CG + BoomerAMG
[]

[Reporters]
  [perf]
    type = PerfGraphReporter
    execute_on = FINAL
  []
[]

[Outputs]
  json = true                       # timing sections land in the JSON reporter output
[]
```

- Transient executioners already work with linear systems (for example
  `test/tests/time_integrators/implicit-euler/ie-linearfv.i`).
- Without a time-derivative kernel, every step reassembles and resolves the same steady system,
  which is exactly the repeated-assembly loop we want.
- Matrix and vector types are selected on the command line: `-mat_type aijkokkos -vec_type kokkos`
  for the device campaign, defaults for the host-PETSc campaign (H8).

### 5.3 Sweep matrix

| Campaign | Variable | Values |
|---|---|---|
| Single-GPU size sweep | cells per GPU | 16k, 64k, 256k, 1M, 4M, 16M, 32M (as memory allows) |
| Mesh type | topology / quality / numbering | hex; tet; skewed hex; scrambled hex |
| Kernel weight | proxy `n_cell_fields` x `flop_factor` | {0, 4, 13, 24} x {0, 50, 200} |
| Strong scaling | GPUs for a fixed 16M-cell mesh | 1, 2, 4, 8 (and 16 if available) |
| Weak scaling | GPUs at 4M cells per GPU | 1, 2, 4, 8, 16 |
| CPU reference | MPI ranks on one full CPU node, CPU objects | same meshes up to node memory |
| Poisson / AMG (H10) | preconditioner policy | rebuild every step; freeze k = 5, 10, 20; GAMG with `reuse_interpolation` |

**Ablation order.** Run the baseline, then each switch alone, then cumulative in the planned order:
fences → slots → `face_gather` → owned COO → memoization → combined hooks → renumbering. This
reveals interactions; for example, removing fences may change which assembly mode wins.

### 5.4 Measurement hygiene

- Use a dedicated GPU node, one MPI rank per GPU, and the same node for all variants of a campaign.
- Lock GPU clocks if allowed (`nvidia-smi -lgc`), or at least record them.
- Discard warm-up steps. The first assembly includes setup, allocation, and COO preallocation;
  report it separately as **setup time**.
- Repeat each configuration in at least 3 separate processes. Report the median and the
  interquartile range of per-assembly time. Treat a difference as real only if it exceeds both
  3% and the IQR.
- Also measure **STREAM bandwidth** on the same GPU (BabelStream with its Kokkos backend). Convert
  face-pass and gather times to achieved bandwidth using a byte model per face and per cell
  (written down per kernel). Compare against miniAero's 34-36% of STREAM as a sanity reference.

### 5.5 Tools

| Question | Tool |
|---|---|
| Time per section | `PerfGraphReporter` JSON (section 5.2), MOOSE perf graph console output |
| Kernel count, gaps between kernels, syncs, memcpy traffic | Nsight Systems (`nsys profile -t cuda,nvtx,mpi`) with the Kokkos Tools NVTX connector |
| Fences per assembly | Kokkos Tools `space-time-stack` |
| Achieved bandwidth, registers, spills, L2 hit rate of face pass and gather | Nsight Compute on representative kernels |
| PETSc time split (`MatSetValuesCOO`, `KSPSolve`, `PCSetUp`, `VecScatter`) | `-log_view` (use `-log_view_gpu_time` only in diagnostic runs; it adds syncs) |
| Device memory high-water mark | Kokkos Tools `memory-hwm` and `nvidia-smi` sampling |

---

## 6. Deliverables and reporting

For each campaign, produce:

1. A table of median per-assembly time per switch combination and size, with the speedup over the
   baseline.
2. Plots:
   - assembly time per cell against cells per GPU (log-log), one line per variant;
   - the H2b crossover heat map over (`n_cell_fields`, `flop_factor`);
   - strong- and weak-scaling efficiency.
3. Kernels and fences per assembly, and H2D/D2H bytes per assembly, per variant.
4. Achieved bandwidth as a fraction of STREAM for the face pass and the gather.
5. A one-line verdict per hypothesis (adopt / reject / inconclusive), applying the rules in
   section 2.

Feed the results back into `kokkos_ns_linearfv_detailed_design.md`: settle R3 (assembly mode),
Q4, Q7, Q8 and Q9, and set the target cells per GPU.

---

## 7. Order of work

| Step | Work | Code needed | Answers |
|---|---|---|---|
| 0 | **Baseline profile of the current code**: the benchmark input above at 3-4 sizes, perf graph, Nsight Systems, `space-time-stack` | Inputs only, plus the finer timing sections (3.2) | Current launch, fence and sync counts; where time goes; whether host PETSc types cause round trips |
| 1 | `sizeof` report; ghost-update microbenchmark | Tiny | H7 / Q7, Q4 / R6 |
| 1b | Mesh dump object + standalone layout mini-app (section 4A) | Medium, independent of MOOSE internals | I3-I7, I9, I10: picks the layout before phase-0 code is written |
| 2 | Per-assembly fence switch | Small (dispatcher) | H1 |
| 3 | Face SoA + `face_atomic` + precomputed slots | Medium (phase-0 infrastructure) | H3; first face-indexed results |
| 4 | `face_gather` + owned COO + memoization + combined hooks | Medium | H2a, H2c, H4, H6 |
| 5 | Proxy flux and gradient objects; renumbering | Medium | H2b, H5 |
| 6 | Multi-GPU scaling; AMG policy study | Inputs only | H9, H10 |

**Step 0 is the most valuable first move.** It needs no design code. It may show that the
dominant cost today is something other than what the design targets, for example the per-assembly
`sync` copies or the host-side PETSc types. That would reorder everything after it.

---

## 8. Environment notes

- **Hardware.** Performance numbers must come from a GPU build (`./configure --with-kokkos` with
  PETSc built with CUDA and Kokkos, see `modules/doc/content/getting_started/installation/install_kokkos.md`)
  on a GPU node.
- **CPU builds.** A `--with-kokkos=cpu` build (OpenMP backend, for example on a laptop) is useful
  only for the correctness gate and determinism checks, not for timing.
- **Thread limits.** The current Kokkos linear FV tests run with `max_threads = 1`. Benchmark runs
  use one MPI rank per GPU and no MOOSE threads.
