#!/usr/bin/env bash
# Regenerates the FVMeshDumper meshes used by the layout mini-app (about 2.1 GB in total), plus
# PETSc matrix/RHS dumps of the small meshes, which the mini-app checks its V0 assembly against.
# Usage: make_meshes.sh [moose_test-opt]   (run inside the moose-next conda env)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
EXE="${1:-$HERE/../../../test/moose_test-opt}"
INPUT="$HERE/../inputs/advdiff_kokkos.i"
cd "$HERE"

dump() {
  local name=$1 ref=$2
  shift 2
  local extra=()
  if [ "$ref" = ref ]; then
    extra=(-ksp_view_mat "binary:$HERE/${name}_A.bin" -ksp_view_rhs "binary:$HERE/${name}_b.bin")
  fi
  "$EXE" -i "$INPUT" steps=1 "$@" UserObjects/dump/type=FVMeshDumper \
    UserObjects/dump/file="$HERE/$name.fvd" UserObjects/dump/variable=u \
    Outputs/file_base="$HERE/tmp_$name" "${extra[@]}" > /dev/null
  echo "wrote $name.fvd"
}

dump hex8_n8 ref n=8
dump hexskew_n8 ref n=8 skew=0.3
dump tet4_n6 ref n=6 elem_type=TET4
dump hex8_n64 noref n=64
dump hex8_n128 noref n=128
dump hexskew_n128 noref n=128 skew=0.3
dump tet4_n40 noref n=40 elem_type=TET4
rm -f tmp_*
