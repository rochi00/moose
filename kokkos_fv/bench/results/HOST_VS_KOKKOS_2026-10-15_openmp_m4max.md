# Host vs Kokkos linear FV assembly: Kokkos OpenMP backend, Apple M4 Max (2026-10-15)

The Kokkos linear FV objects assemble the advection-diffusion-source system 2.5-2.7x faster than
the host linear FV objects on one core, and 2.1-4.7x faster on eight cores (Kokkos threads against
host MPI ranks), with the gap widening with the mesh size. Both assemble the same system:
`bench.py parity` matches them entry by entry.

**Setup.**
- Build: MOOSE `rc/kokkos_moose_navier_stokes` at `14e8322f95` (end of phase 1 step 7), with the
  phase 1 step 8 boundary conditions added but not used by these inputs, `--with-kokkos=cpu`, conda
  env `moose-next`, opt build.
- Machine: Apple M4 Max (10 performance and 4 efficiency cores), 36 GB.
- Problem: `inputs/advdiff_cpu.i` and `inputs/advdiff_kokkos.i`, 3D advection (upwind) - diffusion
  (no nonorthogonal correction) - source on a generated orthogonal hex mesh, Dirichlet and Neumann
  boundary conditions.
- Timings: assembly per step (`LinearSystem::computeLinearSystemInternal`), median of 3 repeats,
  difference method (steps 3 and 13). The linear solve is one Richardson/Jacobi sweep and is not
  part of these numbers.
- Raw data: `sweep_20261015-161519.csv` (Kokkos), `sweep_20261015-162417.csv` (host).

## Assembly time per step

| Cells | Host, 1 rank | Kokkos, 1 thread | Speedup | Host, 8 MPI ranks | Kokkos, 8 threads (1 rank) | Speedup |
|---|---|---|---|---|---|---|
| 32,768 (32^3) | 23.5 ms | 9.37 ms | 2.5x | 6.04 ms | 2.90 ms | 2.1x |
| 262,144 (64^3) | 200 ms | 74.0 ms | 2.7x | 43.6 ms | 16.1 ms | 2.7x |
| 2,097,152 (128^3) | 1557 ms | 602 ms | 2.6x | 479 ms | 102 ms | 4.7x |

Per cell:

| Cells | Host, 1 rank | Kokkos, 1 thread | Host, 8 ranks | Kokkos, 8 threads |
|---|---|---|---|---|
| 32^3 | 718 ns | 286 ns | 184 ns | 88 ns |
| 64^3 | 764 ns | 282 ns | 166 ns | 61 ns |
| 128^3 | 742 ns | 287 ns | 228 ns | 49 ns |

- **One core.** The per-cell cost is flat for both, so the ratio is a property of the code: the
  Kokkos face loop over contiguous arrays against the host's per-face object calls and functor
  evaluations.
- **Eight cores.** Kokkos speeds up 3.2x, 4.6x and 5.9x over its own single thread; the host 3.9x,
  4.6x and 3.2x over its single rank. At 128^3 the host's 8-rank time degrades.
- **Noise.** The 8-core runs at 128^3 vary a lot between repeats: interquartile range 0.29 s on
  the host's 0.48 s median, and 0.05 s on Kokkos' 0.10 s. The 4.7x there is less certain than the
  other ratios. At 32^3 and 64^3 the 8-core interquartile ranges are under 10% of the median for
  the host and 34-39% for Kokkos.

## Memory (maximum per rank)

| Cells | Host, 1 rank | Kokkos, 1 rank | Host, 8 ranks (per rank) | Kokkos, 8 threads |
|---|---|---|---|---|
| 32^3 | 237 MB | 376 MB | 232 MB | 387 MB |
| 64^3 | 1287 MB | 2328 MB | 835 MB | 2260 MB |
| 128^3 | 4805 MB | 7332 MB | 2732 MB | 7003 MB |

On one rank Kokkos uses 1.5-1.8x the host's memory: it keeps the Kokkos mesh, system and FV mesh
copies of the data alongside the MOOSE ones. Eight host ranks together use about 3x the memory of
one Kokkos process with eight threads.

## Not covered

These are CPU-backend numbers. Kokkos on a GPU against the host on CPU cores needs the GPU
campaign; the linear solve, which both hand to the same PETSc solver, is excluded.
