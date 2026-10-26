#!/usr/bin/env bash
# GPU measurement campaign for the open design questions (kokkos_fv/docs/kokkos_ns_linearfv_detailed_design.md,
# section 15.2). Run inside a GPU allocation, from kokkos_fv/bench, after building:
#   - MOOSE (test/moose_test-opt) with --with-kokkos against a GPU PETSc,
#   - the mini-app against the same Kokkos (cmake -S miniapp -B miniapp/build -DKokkos_ROOT=...),
#   - the probe tools/petsc_kokkos_probe.cpp (see its header) as tools/petsc_kokkos_probe,
#   - the event counter tools/kokkos_event_counter.cpp as tools/libkec.so.
#
# Usage: ./run_gpu_campaign.sh [section ...]   (default: all sections except miniapp_gradients,
#                                               in order)
# Sections: correctness determinism miniapp miniapp_gradients overhead phase1 types q4 q8
#           coefficients solvers ranks
#
# Environment:
#   BENCH_EXE      MOOSE executable (default ../../test/moose_test-opt)
#   BENCH_MPIEXEC  MPI launcher accepting -n (default mpiexec; srun works)
#   MESHES         mini-app meshes (default meshes; see meshes/make_meshes.sh)
#   NRANKS         rank counts for the 'q4' and 'ranks' sections (default "1 2 4")
#
# On GPU builds, Kokkos-only linear systems now select -<sys>_mat_type aijkokkos -vec_type kokkos
# by default; host PETSc types are requested explicitly below where they are compared.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
EXE="${BENCH_EXE:-$HERE/../../test/moose_test-opt}"
MPIEXEC="${BENCH_MPIEXEC:-mpiexec}"
export MESHES="${MESHES:-$HERE/meshes}"
NRANKS="${NRANKS:-1 2 4}"
STAMP=$(date +%Y%m%d-%H%M%S)
LOGDIR="$HERE/results/gpu_campaign_$STAMP"
mkdir -p "$LOGDIR"
HOST_TYPES="-vec_type standard -u_sys_mat_type aij"
DEVICE_TYPES="-vec_type kokkos -u_sys_mat_type aijkokkos"
GATHER="Problem/kokkos_linear_fv_assembly=gather"
# Linear FV variables default to the host 'green-gauss' gradient method; the Kokkos one must be
# selected. Command-line blocks are merged into the input, so it needs no input change.
KOKKOS_GRADIENT="FVGradientMethods/kgg/type=KokkosFVGreenGaussGradient Variables/u/gradient_method=kgg"

log() { echo "=== $*" | tee -a "$LOGDIR/campaign.log"; }

# Kokkos assembly against the host assembly, entry by entry: hex, tet and prism, orthogonal and
# skewed, with and without non-orthogonal correction, every advection interpolation, 1 and 3
# ranks (240 cases). On the CPU build every case matches.
correctness() {
  log "Parity of the Kokkos and host assembled systems"
  python3 bench.py parity --interp upwind average minmod vanleer muscl 2>&1 \
    | tee "$LOGDIR/parity.txt"
}

# Run-to-run reproducibility of the assembled system: atomic (default) and gather modes on a
# skewed tet mesh, where atomic assembly differed between runs on the CPU
determinism() {
  log "Determinism, atomic assembly (differing runs expected)"
  python3 bench.py determinism --elem TET4 --n 24 --runs 5 --set skew=0.3 \
    | tee "$LOGDIR/determinism_atomic.txt" || true
  log "Determinism, gather assembly (one distinct system expected)"
  python3 bench.py determinism --elem TET4 --n 24 --runs 5 --set skew=0.3 "$GATHER" \
    | tee "$LOGDIR/determinism_gather.txt" || true
}

# Launch cost (Q14) and the fused face pass (Q13) with the standalone mini-app
miniapp() {
  log "mini-app: launch cost (Q14)"
  miniapp/build/fv_layout_bench --launch-bench | tee "$LOGDIR/launch.csv"
  log "mini-app: assembly layouts incl. fused V2F (Q13)"
  SCALING_THREADS="" miniapp/run_layout_bench.sh 20 | tee -a "$LOGDIR/campaign.log"
}

# Gradient layouts (Q11, Q12, Q18), answered on an A100; rerun only on a different GPU
miniapp_gradients() {
  log "mini-app: gradient layouts G0-G3, G0C, G1C, G1CS (Q11, Q12)"
  THREADS=1 miniapp/run_gradient_bench.sh 20 3 | tee -a "$LOGDIR/campaign.log"
}

# Per-assembly host overhead after the fence/allocation changes (H1), and per-section times
overhead() {
  log "MOOSE steady sweep, device types (default on GPU) and host types"
  python3 bench.py sweep --impl kokkos --n 16 32 64 128 --repeats 3 --label gpu_device_types
  python3 bench.py sweep --impl kokkos --n 16 32 64 128 --repeats 3 --label gpu_host_types \
    --petsc-args "$HOST_TYPES"
  log "Kokkos event counts per assembly (difference of 2 and 12 assemblies)"
  for s in 2 12; do
    KOKKOS_TOOLS_LIBS="$HERE/tools/libkec.so" KEC_OUT="$LOGDIR/events_$s.txt" \
      "$EXE" -i inputs/advdiff_kokkos.i n=8 steps=$s Outputs/file_base="$LOGDIR/events_$s" > /dev/null
  done
  log "Transient per-step times: device types + BiCGStab, and host types + GMRES"
  python3 bench.py transient --variants cpu kokkos --n 64 --steps 20 --repeats 3 \
    --set "Executioner/petsc_options_value=bcgs jacobi" --label gpu_device_bcgs
  python3 bench.py transient --variants kokkos --n 64 --steps 20 --repeats 3 \
    --petsc-args "$HOST_TYPES" --label gpu_host_gmres
  log "Steady sweep, gather assembly"
  python3 bench.py sweep --impl kokkos --n 64 128 --repeats 3 --label gpu_gather --set "$GATHER"
}

# Phase 1 numerics on a skewed hex mesh: gradients (column cpu_gradients, which times every
# gradient method, host or Kokkos), non-orthogonal correction, limited advection and a property
# from a KokkosParsedFunctorMaterial. Minmod and Van Leer always request the host 'green-gauss'
# gradient; MUSCL takes a gradient method, here the Kokkos one with the Venkatakrishnan limiter.
phase1() {
  local skewed="skew=0.3"
  local sweep="python3 bench.py sweep --impl kokkos --n 64 128 --repeats 3"
  log "Phase 1: skewed baseline (no gradients)"
  $sweep --label p1_base --set $skewed
  log "Phase 1: non-orthogonal correction, host and Kokkos gradients"
  $sweep --label p1_noc_host_grad --set $skewed noc=true
  $sweep --label p1_noc_kokkos_grad --set $skewed noc=true $KOKKOS_GRADIENT
  log "Phase 1: Van Leer (host gradient) and MUSCL (limited Kokkos gradient)"
  $sweep --label p1_vanleer_host_grad --set $skewed interp=vanleer
  $sweep --label p1_muscl_kokkos_grad --set $skewed interp=muscl \
    FVGradientMethods/kggv/type=KokkosFVGreenGaussGradient \
    FVGradientMethods/kggv/limiter=venkatakrishnan FVInterpolationMethods/muscl/gradient_method=kggv
  log "Phase 1: diffusion coefficient from a parsed functor material of u"
  $sweep --label p1_functor_material --set $skewed \
    FunctorMaterials/kmat/type=KokkosParsedFunctorMaterial \
    "FunctorMaterials/kmat/expression=1+0.5*x+0.1*u" FunctorMaterials/kmat/functor_names=u \
    FunctorMaterials/kmat/property_name=k_mat LinearFVKernels/diffusion/diffusion_coeff=k_mat
}

# The default device-type selection (section 5 of findings_2026-10-14.md)
types() {
  log "Device PETSc types selected by default"
  "$EXE" -i inputs/advdiff_kokkos.i n=8 steps=1 Outputs/console=true \
    Outputs/file_base="$LOGDIR/types" -ksp_view 2>&1 | grep -E "keeping its matrix|type: (seq|mpi)aij" \
    | sort | uniq -c | tee "$LOGDIR/types.txt"
}

# Q4: are ghost updates of libMesh ghosted vectors with Kokkos vector types done on the device?
q4() {
  for np in $NRANKS; do
    [ "$np" -lt 2 ] && continue
    log "Q4: -log_view with $np ranks, device types"
    $MPIEXEC -n "$np" "$EXE" -i inputs/advdiff_transient_kokkos.i n=64 steps=5 \
      Outputs/file_base="$LOGDIR/q4_np$np" $DEVICE_TYPES -log_view \
      > "$LOGDIR/q4_np$np.log" 2>&1
    # GPU-host copy counts per event; VecScatter/VecGhostUpdate rows show where ghosts move
    grep -E "^(Event |VecScatter|SFBcast|SFReduce|SFPack|SFUnpack|VecGhost|MatSetValuesCOO|VecSetValues)" \
      "$LOGDIR/q4_np$np.log" | tee -a "$LOGDIR/q4_summary.txt" || true
  done
}

# Q8: does PETSc run on Kokkos's default execution space instance?
q8() {
  log "Q8: PETSc/Kokkos execution space probe"
  tools/petsc_kokkos_probe | tee "$LOGDIR/q8.txt"
}

# Q15: cost of the parsed-function interpreter on device (coefficient pass with a parsed
# coefficient against a constant one)
coefficients() {
  log "Q15: parsed diffusion coefficient vs constant"
  python3 bench.py sweep --impl kokkos --n 64 128 --repeats 3 --label gpu_k_parsed
  python3 bench.py sweep --impl kokkos --n 64 128 --repeats 3 --label gpu_k_constant \
    --set "Functions/k_fn/expression=1.5"
}

# Q16: linear solvers on device. Poisson (no advection: symmetric, CG + AMG) and the
# advection-diffusion problem (nonsymmetric)
solvers() {
  local poisson="LinearFVKernels/advection/velocity=0 0 0"
  log "Q16: Poisson, CG with Jacobi, hypre BoomerAMG and GAMG"
  python3 bench.py transient --variants kokkos --n 64 128 --steps 10 --repeats 3 \
    --set "$poisson" "Executioner/petsc_options_value=cg jacobi" --label poisson_cg_jacobi
  python3 bench.py transient --variants kokkos --n 64 128 --steps 10 --repeats 3 \
    --set "$poisson" "Executioner/petsc_options_iname=-ksp_type -pc_type -pc_hypre_type" \
    "Executioner/petsc_options_value=cg hypre boomeramg" --label poisson_cg_boomeramg
  python3 bench.py transient --variants kokkos --n 64 128 --steps 10 --repeats 3 \
    --set "$poisson" "Executioner/petsc_options_iname=-ksp_type -pc_type -pc_gamg_reuse_interpolation" \
    "Executioner/petsc_options_value=cg gamg true" --label poisson_cg_gamg_reuse
  log "Q16: advection-diffusion, Krylov methods with fewer reductions"
  for ksp in gmres bcgs pipebcgs; do
    python3 bench.py transient --variants kokkos --n 64 128 --steps 10 --repeats 3 \
      --set "Executioner/petsc_options_value=$ksp jacobi" --label advdiff_${ksp}_jacobi
  done
}

# Q17: ranks per GPU and multi-GPU scaling (set NRANKS to match the node's GPU count, and once
# with more ranks than GPUs to measure oversubscription)
ranks() {
  log "Q17: rank scaling, device types + BiCGStab"
  python3 bench.py transient --variants kokkos --n 128 --np $NRANKS --steps 10 --repeats 3 \
    --set "Executioner/petsc_options_value=bcgs jacobi" --label ranks_bcgs
}

sections=("$@")
[ ${#sections[@]} -eq 0 ] &&
  sections=(correctness determinism miniapp overhead phase1 types q4 q8 coefficients solvers ranks)
for s in "${sections[@]}"; do
  "$s"
done
log "Logs in $LOGDIR; CSVs in results/"
