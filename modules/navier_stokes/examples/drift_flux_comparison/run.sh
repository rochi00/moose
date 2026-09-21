#!/bin/sh
# Run one of the comparison inputs and record its wall time, step count, failed solves and, for
# the pipe, the drift-flux parameters at the z/D = 54 station. Serial by default; set NP for MPI.
#
#   ./run.sh <input.i> [extra MOOSE command line arguments]
#
#   NP=4 ./run.sh pipe/hibiki-50mm-pipe.i Mesh/uniform_refine=1
#   NP=4 ./run.sh pipe/newton-pipe.i MUMPS
#   ./run.sh tjunction/pimple.i
#
# The word MUMPS is replaced by the options that hand the Newton factorization to MUMPS, which
# is what to use on more than one rank. The executable is looked for next to this directory
# (the usual in-tree build) and then on the PATH; override with NAVIER_STOKES_EXE. MPI is run
# through the mpiexec on the PATH; override with MPIEXEC.
set -e
here=$(cd "$(dirname "$0")" && pwd)
input=$1; shift
[ -n "$input" ] || { sed -n '2,15p' "$0"; exit 1; }

exe=${NAVIER_STOKES_EXE:-$here/../../navier_stokes-opt}
[ -x "$exe" ] || exe=$(command -v navier_stokes-opt) || { echo "navier_stokes-opt not found"; exit 1; }
mpi=${MPIEXEC:-mpiexec}
np=${NP:-1}

args=""
for a in "$@"; do
  if [ "$a" = MUMPS ]; then
    args="$args Preconditioning/SMP/petsc_options_iname='-pc_type -pc_factor_mat_solver_type -pc_factor_shift_type'"
    args="$args Preconditioning/SMP/petsc_options_value='lu mumps NONZERO'"
  else
    args="$args $a"
  fi
done

tag=$(basename "$input" .i)_np$np
cd "$(dirname "$input")"
start=$(date +%s)
eval $mpi -n $np "$exe" -i "$(basename "$input")" $args Outputs/file_base=$tag Outputs/perf_graph=true > $tag.log 2>&1
rc=$?
wall=$(( $(date +%s) - start ))
station=$(tail -1 $tag.csv 2>/dev/null | awk -F, 'NF > 20 {printf "C0_54=%.4f Vgj_54=%.4f", $3, $6}')
echo "$tag rc=$rc wall=${wall}s steps=$(grep -c 'Time Step' $tag.log) failed=$(grep -c 'Solve Did NOT' $tag.log) $station" | tee -a "$here/results.log"
