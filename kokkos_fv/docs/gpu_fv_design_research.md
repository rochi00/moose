# Scalable GPU CFD Designs: Lessons for Kokkos Navier-Stokes Linear FV

This is a survey of production and research GPU CFD codes, done to inform
`kokkos_ns_linearfv_design.md`. Pitfalls are referenced by their numbers from the earlier risk
review: P1 (launch/fence latency), P2 (halo count), P3 (reductions), P4 (pressure AMG),
P5 (small sequential solves), P6/P7 (gather/scatter locality), P8 (numbering), P9 (registers),
P10 (functor size), P12 (memory), P14 (ghost correctness), P15 (load balance), P17 (outer
iteration count).

**Confidence tags.** Claims were gathered by web research in October 2026. Untagged claims come
from sources that were read directly. **[S]** marks claims seen only in a search snippet or
abstract. **[U]** marks claims that are uncertain or inferred. Verify any number before quoting it
externally.

---

## 1. Key takeaways

1. **Once assembly is on the GPU, linear solves and launch/sync overhead dominate. Assembly
   does not.**
   - SPUMA (a full OpenFOAM GPU port) measured, per SIMPLE iteration at 65M cells on 4 GPUs:
     momentum (U) solve 0.55 s, pressure solve 0.30 s, pressure assembly 0.04 s.
     https://arxiv.org/html/2512.22215v1
   - This corrects our earlier assumption that "pressure is 50-80% of time". Measured pressure
     shares are about 5-28% (SPUMA) and 30-40% (Fluent segregated,
     https://www.nvidia.com/en-au/data-center/gpu-accelerated-applications/ansys-fluent/).
     **Momentum solves can cost more than pressure.**
2. **Synchronizing after every kernel is the main avoidable cost.**
   - In SPUMA, `cudaDeviceSynchronize` after every kernel (about 53k calls in 5 iterations) took
     73-92% of CUDA API time.
   - Kokkos-MOOSE fences after every dispatch, which is the same pattern (P1).
   - A synchronous launch costs about 7-13 µs, against about 2-5 µs asynchronous.
     https://arxiv.org/html/2406.08646v2
3. **Strong scaling ends at about 1-10M cells per GPU no matter how the code is designed.**
   - SPUMA: about 65% efficiency at about 8M cells/GPU.
   - OGL: utilization is good only above about 1M DOF per GPU.
   - NekRS: 80% efficiency at about 2M points per V100.
   - Commercial production runs use 8-30M cells per GPU (STAR-CCM+ 20-27M, Fluent 8-31M).
   - Design for at least about 2-5M cells per GPU, and treat memory capacity as the second
     scaling axis (P12).
4. **No single way of avoiding write races wins everywhere. Benchmark per kernel, and keep the
   choice switchable.**
   - Atomics: Nalu-Wind, Albany, SPUMA, code_saturne (CUDA), and Zhang et al. (FGCS 2023), where
     atomics beat coloring at every size on V100.
   - Two-step gather (store face values, then gather per cell, which is our design): beat atomics
     in miniAero; AMD moved OpenFOAM in this direction in 2026 ("face-to-cell transformations").
   - Coloring: SPARC. It lost to atomics wherever it was measured directly.
   - Choosing the loop mode per kernel gave 15-20% (Weng 2021 [S]).
   - Determinism is a separate reason to prefer gathers. One production code added cell coloring
     for bitwise reproducibility (Appl. Sci. 2023 [S]).
5. **Reuse is what makes AMG on GPU affordable.**
   - One GPU AMG setup costs about as much as one solve.
   - A device memory pool alone cut hypre setup 4.6x.
   - Nalu-Wind rebuilds its preconditioner only every 100 solves in its GPU configuration.
   - Nalu-Wind solves u, v and w with **one** matrix and preconditioner (P4, P5).
6. **Everything in the iteration loop must stay on the device.**
   - Codes that offloaded only the solver: "no or only limited advantages" (EXASIM report).
   - Incomplete ports (boundary conditions, turbulence) gave 0.4x in the stdpar study [S].
   - Page migration took more than 65% of runtime on discrete GPUs with unified memory (AMD).
   - Any CPU object left in a Kokkos FV system's loop wipes out the gain.

---

## 2. Codes surveyed

### Kokkos-based

| Code | Mesh / data | Assembly | Solver | Lessons |
|---|---|---|---|---|
| **Nalu-Wind / ExaWind** (renamed Kynema-UGF 2026; node-centered edge-based FV, incompressible, segregated) | STK NGP buckets: contiguous same-topology chunks, one Kokkos team per bucket | Atomics into a device CSR value array. The sparsity pattern is built once. Buffers are laid out as `[owned nnz \| shared nnz]`; shared entries go through hypre `AddToValues` and are merged on device by sort+reduce | hypre GMRES. Pressure: BoomerAMG with PMIS, aggressive coarsening, MM-ext interpolation, two-stage GS or l1-Jacobi. Momentum/scalars: single-level smoother. `recompute_preconditioner_frequency: 100`. `HypreUVWSolver` solves u, v, w with one matrix and preconditioner | Scaled to more than 4000 V100s. AMG setup and hypre assembly were the bottlenecks. Strong-scaling slope -0.47 on GPU vs -0.98 on CPU. The MPI stack mattered: Eagle beat Summit by 15% with half the GPUs. https://docs.nlr.gov/docs/fy22osti/81212.pdf , https://github.com/kynema/kynema-ugf , https://github.com/Exawind/nalu-wind (`reg_tests/hypre_settings/hypre_blade_resolved.yaml`) |
| **SPARC** (Sandia, cell-centered compressible FV) | One 2D View `[entities x vars]` per block, LayoutLeft so entities vary fastest (coalesced) | 3D graph coloring, kept until the mesh changes; entities reordered within each color; "atomics found slower" (2018, older GPUs) | Line/block Jacobi (Kokkos Kernels) | Small boundary kernels "do not scale" because of launch time. End-to-end time was similar across architectures. Long builds from template composition. https://www.iccfd.org/iccfd10/papers/ICCFD10-079-Paper.pdf |
| **Albany / MALI** (Kokkos, Phalanx, Tpetra) | Owned+shared rows in one `FECrsMatrix` | Atomic scatter, then export | Tpetra/Ifpack2 | One matrix with owned and shared rows: 2x less memory and 2.1x faster assembly. **Boundary data aligned 1:1 with the side list: 8.7x on GPU** (it was index-indirected before). https://arxiv.org/pdf/2204.04321 |
| **miniAero** (Kokkos unstructured FV) | Face loop | Atomics or two-step gather, selected at compile time. Shared-memory hierarchical coloring with partition reordering beat both (by 92% over two-step gather, 158% over atomics on V100) | n/a | The best variant still reached only 34-36% of STREAM bandwidth. https://arxiv.org/pdf/1802.03749 |
| **Uintah** | Structured AMR | Kokkos loops inside an asynchronous task runtime | n/a | Bulk-synchronous Kokkos loops replaced by asynchronous tasks on streams to hide latency [S]. https://www.sci.utah.edu/publications/Hol2019a/UUSCI-2019-001.pdf |

### GPU ports of segregated finite volume solvers

| Code | Design | Lessons |
|---|---|---|
| **SPUMA** (OpenFOAM, CUDA/HIP, 2025) | Native LDU format kept; face loops with atomics; managed memory with a pool; CRTP lambda executor | Memory pool: about 4-6x. Synchronization after every kernel dominated API time. Gradient kernels were among the slowest because they are atomic-heavy. Native GAMG built too many tiny levels. AmgX rebuilt its hierarchy every step. https://arxiv.org/html/2512.22215v1 |
| **AMD OpenFOAM offload** (OpenMP, unified memory) | Loops offloaded only above a size cutoff; Umpire pools | MI300A was 4x H100 because discrete GPUs lost more than 65% of time to page migration. The 2026 follow-up names launch overhead, atomics and memory management as the bottlenecks, with fixes: kernel fusion via expression templates, face-to-cell loops, pools. https://arxiv.org/pdf/2405.00436 |
| **code_saturne** (EDF) | `cs_dispatch` lambdas with `parallel_for_i_faces` / `_b_faces`; `cs_dispatch_sum` is SIMPLE (race-free face numbering on CPU) or ATOMIC (CUDA); one kernel body for both | Renumbering options: Morton, Hilbert, METIS. An earlier version fused kernels with template pack expansion to cut launches (SC18 [S]). https://raw.githubusercontent.com/code-saturne/code_saturne/master/src/base/cs_dispatch.h |
| **OGL / Ginkgo, amgx4Foam, PETSc4FOAM** (solver-only offload) | LDU→CSR/COO **permutation computed once**; only values updated each step | Doing this conversion naively dominated runtime (MetaX). GPU-aware MPI: 25-50%. Oversubscribing a GPU without repartitioning: 0.007x. https://arxiv.org/html/2510.08536v1 , https://doi.org/10.1007/s11012-024-01806-1 |
| **Fluent native GPU, STAR-CCM+ GPU** | Built from scratch so everything stays on the GPU, including monitors | One A100 ≈ 540 cores at 8-31M cells/GPU (Fluent). STAR-CCM+ production runs at 20-27M cells/GPU. Memory, not compute, sets the GPU count. |
| **Zhang et al. FGCS 2023** (MPI+CUDA unstructured FV) | Atomics; renumbered cells and faces; halo pack/unpack on the GPU | Renumbering: 1.4-1.5x per kernel. Strong scaling 42% at 200 GPUs; weak scaling 82%. https://xianweiz.github.io/doc/papers/cfd_fgcs23.pdf |
| **ICSFoam + AmgX** (block-coupled p-U) | Coupled pressure-velocity solve | More than 4x end-to-end over CPU SIMPLE on 45M cells; linear-solver share fell from 80% to 20-25%. Almost no gain at 128M cells. https://arxiv.org/html/2403.07882v1 |

---

## 3. Techniques, mapped to our pitfalls

### P1. Launch latency and fences

| Technique | Evidence | Applicability |
|---|---|---|
| **Remove per-dispatch fences.** Fence only before host reads or non-stream-ordered MPI; use instance fences | SPUMA (73-92% of API time); a sync costs 2.5-4x a launch | **Highest priority.** Needs a "no fence" option in `Dispatcher` for FV passes. |
| **Fuse kernels**: one face kernel computing convection + diffusion + non-orthogonal correction, all boundary kinds in one pass | AMD 2026, code_saturne SC18, SPARC boundary kernels | Fuse within one object's hooks first. Fusing across objects needs a composition mechanism in the dispatcher. |
| **Kokkos Graph / CUDA graphs** | Break-even at about 40 replays; 1.38x at 500 replays on tiny kernels; about 60 ns saved per node; **MPI and PETSc calls cannot be inside a graph** [U] | Last resort. Only for fixed kernel sequences between MPI points. https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Graph.html |
| **Execution-space instances** (`partition_space`) for overlap | Kokkos docs; Uintah | Useful for overlapping interior work with the halo exchange (P2). |

### P2. Halo exchanges

| Technique | Evidence | Applicability |
|---|---|---|
| **Pack fields into one ghosted block vector** (`VecCreateGhostBlock`, bs = number of fields) | PetscSF sends multi-field units as one message. Each separate exchange costs at least 35-45 µs (Summit) | Store gradient components (and HbyA, Ainv) as one blocked vector: 20+ exchanges become about 3-5. https://petsc.org/main/manualpages/Vec/VecCreateGhostBlock/ |
| **Overlap the interior pass with the exchange** | NekRS: 10-15% saved | Matches our face-range order (processor-boundary faces last). |
| **GPU-aware MPI with device buffers** | +25-50% (OGL), but it was slower or buggy on Summit and Perlmutter (PETSc paper) | Make it configurable and measure per machine. |
| **Two ghost layers** to compute ghost gradients locally | Used in a multi-GPU gas-kinetic FV code; no unstructured measurements [U] | An alternative to exchanging gradients. Must respect P14: ghost cells need their complete face set. |
| **kokkos-comm** | "Not ready for general use" | Do not depend on it. |

**PETSc ghosted vectors on device (Q4).** The PETSc source shows that `VECMPIKOKKOS` +
`VecMPISetGhost` builds a ghosted vector whose local form is device-resident, with ghost updates
through PetscSF. However, a February 2026 petsc-users reply says device ghost vectors are not
supported [U, conflicting]. Plain `VecCreateGhost` (used by libMesh) builds a host `VECMPI`.
**Q4 is still open; settle it with a small Nsight test.**

### P3. Global reductions

| Technique | Evidence | Applicability |
|---|---|---|
| Fixed iteration counts or lagged norms for momentum (`-ksp_norm_type none`); compute SIMPLE residual norms less often | PETSc guidance | Conflicts with the executioner's convergence test, which uses each system's initial residual. The norm-free variant would need to be optional. |
| Pipelined Krylov (`pipecg`, `pgmres`) | About 2x at scale on CPU; no GPU data; needs asynchronous MPI progress [U] | Experiment later. |
| Device-resident scalars (`VecDotAsync`) | Up to 80% for small CG; KSP integration incomplete | Watch, but do not plan on it. |

### P4/P5. Linear solvers

- **Pressure.**
  - On a device matrix, PETSc already defaults hypre to GPU-friendly settings (from `ihypre.c`):
    PMIS coarsening, ext+i interpolation, l1-Jacobi smoothing, keep-transpose.
  - Add aggressive coarsening (`agg_nl 2`), `P_max 3`, `truncfactor 0.1`, and a strong threshold
    swept over 0.25-0.6. MOOSE's 3D default of 0.7 is tuned for CPU.
  - Use **hypre ≥ 2.31**. Older versions re-convert the matrix on every setup.
  - Enable a **device memory pool** (`HYPRE_SetGPUMemoryPoolSize` or Umpire).
  - Reuse:
    - Freeze the preconditioner (`KSPSetReusePreconditioner`) and rebuild every k outer
      iterations or when the iteration count grows (Nalu-Wind uses k = 100).
    - Alternatively, use **GAMG with `-pc_gamg_reuse_interpolation`**, which recomputes only the
      Galerkin products on device. That is a true coefficient-only resetup, which hypre does not
      offer.
  - Benchmark both.
- **Momentum, energy, k/ε.**
  - GMRES or BiCGStab with Jacobi or l1-Jacobi (Nalu-Wind: a single-level BoomerAMG smoother).
  - No triangular-solve ILU. If ILU is needed, use hypre ILU with iterative triangular solves,
    which was 5x faster.
  - Loose relative tolerance (1e-2 to 1e-3) on non-final outer iterations.
- **One matrix and one preconditioner for u, v, w** (Nalu-Wind `HypreUVWSolver`). This amortizes
  setup 3x and is the low-risk version of the momentum-sharing idea. Avoid `KSPMatSolve` on GPU
  for now (NaNs and extra copies reported).
- Recommended option sets are in section 5.

### P6/P7/P8. Locality

- **Renumber cells (RCM, Hilbert or Morton) and sort faces by owner cell.**
  - Gains measured: 1.1-1.75x in OP2, 1.4-1.5x per kernel in Zhang.
  - Hilbert showed no advantage over Morton in an LBM study.
  - Renumbering also improves AMG and ILU quality.
- **SoA with the entity index varying fastest** (SPARC's LayoutLeft for `(entity, var)`). This
  matches `Moose::Kokkos::Array`'s default.
- **Boundary data aligned 1:1 with the boundary-face list** (Albany: 8.7x on GPU). Our contiguous
  boundary ranges with `bc_face_index = f - range.begin` already do this.
- **Precomputed COO/CSR slots** instead of runtime column search. Nalu-Wind still searches a
  sorted row; the FV stencil is fixed, so we can precompute slots.

### P9/P10. Registers and functor size

- **Functor size.**
  - In current Kokkos CUDA with `__grid_constant__` (Volta+, CUDA 12.1+), functors up to 32 KB
    are passed as kernel parameters.
  - Without it, functors from 512 B to 32 KB use **constant memory**. That path takes a mutex, an
    event sync and a staging copy on **every launch**, a host-blocking serialization point.
  - **Action:** check which path the Kokkos that PETSc downloads (`4.7.04`) uses for MOOSE's
    toolchain. Measure `sizeof` of Kokkos FV objects, and keep device payloads small (POD
    structs of arrays, not whole MOOSE objects).
  - Source: https://raw.githubusercontent.com/kokkos/kokkos/develop/core/src/Cuda/Kokkos_Cuda_KernelLaunch.hpp
- **Registers.** `Kokkos::LaunchBounds` only helps when set tighter than the compiler's own
  choice. Check with `-Xptxas -v` or Nsight Compute, and split heavy face kernels where they spill.

### P12. Memory

- **Persistent allocation and pools.** Pools gave SPUMA 4-6x and hypre setup 4.6x. Audit the
  Kokkos path for per-iteration allocations (for example `Vector::create` during host-to-device
  sync when PETSc vectors live on the host [U]).
- **Avoid managed or unified memory on discrete GPUs.**
- **One copy of each matrix.** Owned rows only (our D6), and free the COO index arrays after
  preallocation. Avoid keeping a Kokkos CSR value array alongside PETSc's own (Albany's 2x
  saving came from removing a duplicate matrix).

### P15. Ranks per GPU

- Run one MPI rank per GPU. Oversubscribing without aggregation was catastrophic (OGL 0.007x).
- Kokkos-MOOSE's device selection assumes several ranks per GPU is normal. For the FV path,
  document one rank per GPU, or rely on MPS.

### P17. Algorithmic ceiling

- Block-coupled p-U solvers gave more than 4x end-to-end on GPU over CPU SIMPLE (ICSFoam), and
  Fluent's coupled solver benefits most from GPUs.
- Keep the face SoA and COO slot design general enough to assemble a 4x4 block-coupled system
  later.
- Mixed precision: low value (1.04-1.1x for AMG); skip it.

---

## 4. Proposed revisions to `kokkos_ns_linearfv_design.md`

| # | Revision | Driven by |
|---|---|---|
| R1 | Add a **no-fence dispatch** option and a single fence per linear-system assembly; make this phase 0 work | SPUMA, P1 |
| R2 | The face pass writes into **per-element-side scratch** `(side, cell)`. The gather writes each cell's whole row **in COO slot order**, so `MatSetValuesCOO` gets a device array laid out exactly like the preallocated pattern. There are no scattered writes and no MPI phase (owned rows only) | miniAero, PETSc COO paper, P6/P7 |
| R3 | Keep an **atomic face-scatter variant** behind a switch, and choose per kernel after benchmarking | Zhang, Nalu-Wind, Weng |
| R4 | **Renumber cells** (RCM or Morton) when the Kokkos mesh is built, and sort faces by owner within each range | OP2, Zhang, P8 |
| R5 | Fuse all of a variable's **boundary-condition passes** into one kernel over all its boundary ranges | SPARC, P1 |
| R6 | Store gradients, `HbyA` and `Ainv` as **blocked ghosted vectors**, one exchange per stage. Overlap the interior face pass with the exchange | PetscSF, NekRS, P2 |
| R7 | **One momentum matrix and one preconditioner** for u, v, w (separate right-hand sides), modelled on Nalu-Wind's UVW solver | Nalu-Wind, P5 |
| R8 | Pressure: hypre ≥ 2.31 with GPU defaults plus aggressive coarsening, a memory pool, and preconditioner freezing every k iterations. Benchmark against GAMG with `reuse_interpolation` | hypre and ECP data, Nalu-Wind, P4 |
| R9 | **Error, or at least warn, if any CPU object participates** in a Kokkos FV system's assembly or executes on its per-assembly `execute_on` | EXASIM, stdpar, AMD |
| R10 | Measure **functor size** and the launch path in phase 0; keep device payloads small | Kokkos launch source, P10 |
| R11 | Correct section 8: the pressure share is about 5-40% once assembly is ported, momentum solves can dominate, and the target is at least 2-5M cells per GPU | SPUMA, Fluent, NekRS |
| R12 | Add NVTX / Kokkos Tools regions per SIMPLE stage from phase 0, and profile before tuning | EDF, SPUMA |

---

## 5. Starting PETSc option sets (untested in MOOSE)

Option names were checked against PETSc `ihypre.c` by the research agent. Numeric values are
starting points to tune, not recommendations.

**Pressure, hypre:**
```
-ksp_type cg -ksp_norm_type unpreconditioned
-pc_type hypre -pc_hypre_type boomeramg
-pc_hypre_boomeramg_coarsen_type PMIS -pc_hypre_boomeramg_interp_type ext+i
-pc_hypre_boomeramg_relax_type_all l1scaled-Jacobi -pc_hypre_boomeramg_no_CF
-pc_hypre_boomeramg_agg_nl 2 -pc_hypre_boomeramg_agg_num_paths 1
-pc_hypre_boomeramg_P_max 3 -pc_hypre_boomeramg_truncfactor 0.1
-pc_hypre_boomeramg_strong_threshold 0.25
-pc_hypre_boomeramg_keeptranspose 1
```

**Pressure, GAMG with coefficient-only resetup:**
```
-pc_type gamg -pc_gamg_type agg -pc_gamg_reuse_interpolation true
-pc_gamg_aggressive_coarsening 1 -pc_gamg_threshold 0.02
-pc_gamg_process_eq_limit 50000 -pc_gamg_coarse_eq_limit 2000
-mg_levels_ksp_type chebyshev -mg_levels_pc_type jacobi -mg_levels_ksp_max_it 2
```

**Momentum and scalars:**
```
-ksp_type bcgs -pc_type jacobi        (or gmres with -ksp_gmres_restart 30)
```

All of these assume `-mat_type aijkokkos -vec_type kokkos`, so matrices and vectors stay on the
device.
