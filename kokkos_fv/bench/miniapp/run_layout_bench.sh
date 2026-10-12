#!/usr/bin/env bash
# Layout benchmark sweep for fv_layout_bench. Writes one CSV with a row per
# (mesh, threads, cell order, face order, variant).
#   usage: ./run_layout_bench.sh [reps]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/build/fv_layout_bench"
MESHES="$HERE/../meshes"
OUT="$HERE/../results/layout_$(date +%Y%m%d-%H%M%S).csv"
REPS="${1:-10}"
export OMP_PROC_BIND=spread OMP_PLACES=threads

echo "mesh,threads,cell_order,face_order,variant,max_rel_diff_vs_V0,distinct_hashes,face_ms,boundary_ms,cell_ms,total_ms,ns_per_cell" > "$OUT"

run() {
  local mesh=$1 threads=$2 co=$3 fo=$4
  "$BIN" --mesh "$MESHES/$mesh.fvd" --cell-order "$co" --face-order "$fo" --reps "$REPS" \
    --kokkos-num-threads="$threads" 2>/dev/null |
    grep -E '^V[0-9]' | sed "s/^/$mesh,$threads,$co,$fo,/" >> "$OUT"
  echo "done $mesh t=$threads $co/$fo"
}

for mesh in hex8_n64 hex8_n128 hexskew_n128 tet4_n40; do
  # full ordering factorial at 8 threads
  for co in native morton rcm scrambled; do
    for fo in native owner; do
      run "$mesh" 8 "$co" "$fo"
    done
  done
  # thread scaling for the baseline ordering and the proposed ordering
  for threads in 1 4; do
    run "$mesh" "$threads" native native
    run "$mesh" "$threads" morton owner
  done
done
echo "wrote $OUT"
