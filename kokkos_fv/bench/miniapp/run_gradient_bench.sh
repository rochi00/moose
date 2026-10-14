#!/usr/bin/env bash
# Green-Gauss gradient sweep for fv_layout_bench. Writes one CSV with a row per
# (mesh, threads, cell order, face order, process repeat, variant).
#   usage: ./run_gradient_bench.sh [reps] [process_repeats]
# Meshes are read from $MESHES (default ../meshes; see ../meshes/make_meshes.sh).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/build/fv_layout_bench"
MESHES="${MESHES:-$HERE/../meshes}"
OUT="$HERE/../results/gradient_$(date +%Y%m%d-%H%M%S).csv"
REPS="${1:-20}"
REPEATS="${2:-3}"
export OMP_PROC_BIND=spread OMP_PLACES=threads

echo "mesh,cells,threads,cell_order,face_order,repeat,variant,max_rel_diff_vs_ref,distinct_hashes,ms,ns_per_cell" > "$OUT"

run() {
  local mesh=$1 threads=$2 co=$3 fo=$4 rep=$5
  local out cells
  out=$("$BIN" --mesh "$MESHES/$mesh.fvd" --cell-order "$co" --face-order "$fo" --reps "$REPS" \
    --variants none --grad-variants G0,G0C,G1,G1C,G1CS,G2,G3 --kokkos-num-threads="$threads")
  cells=$(echo "$out" | sed -n 's/.*cells=\([0-9]*\).*/\1/p' | head -1)
  echo "$out" | grep -E '^G[0-9]' | sed "s/^/$mesh,$cells,$threads,$co,$fo,$rep,/" >> "$OUT"
  echo "done $mesh t=$threads $co/$fo rep=$rep"
}

for mesh in hex8_n64 hex8_n128 hexskew_n128 tet4_n40; do
  # THREADS applies to host backends only; use THREADS=1 on a GPU
  for threads in ${THREADS:-1 4 8}; do
    for order in "native native" "rcm owner" "morton owner" "scrambled native"; do
      for rep in $(seq 1 "$REPEATS"); do
        run "$mesh" "$threads" $order "$rep"
      done
    done
  done
done
echo "Wrote $OUT"
