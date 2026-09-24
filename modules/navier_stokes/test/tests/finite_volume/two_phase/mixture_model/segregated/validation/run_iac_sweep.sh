#!/bin/bash
# Regenerates the interfacial area validation set of the Hibiki interfacial area comparison.
#
# The v5 prefix marks the station definitions of this input: the outlet average is taken at the
# measured station z/D = 53.5, not at 54. The v4 results were read half a diameter downstream of
# the probe, a bias of about 0.1 per cent that always favoured the model.
# Three dispersed phase density models x four void fractions.
#
# 2000 steps at dt = 0.02 is 40 s, about seven residence times. The outlet station settles at about
# 4.4, far later than the inlet ones, so this carries a margin of better than 1.5x; the earlier
# 4000-step sweep was more than was needed.
#
# The relative linear tolerances are set explicitly. Left at their 1e-5 default the converged answer
# is only pinned to about 1e-4, because the linear solver stops iterating once the relative residual
# falls below it and the outer sweep then cannot improve further. Tightening the absolute tolerances
# alone does nothing: the relative criterion binds first.
EXE=/Users/rochishnuchowdhury/Work/INL/moose_two_phase/moose/modules/navier_stokes/navier_stokes-opt
IN=hibiki-iac-compressible.i
# The swarm exponent is 0.25, not the 0.75 of Ishii's correlation itself. The closure forms its
# relaxation time from the continuous phase viscosity, so the factor of (1 - alpha) that the
# mixture density puts into the buoyancy factor survives into the Stokes-limit slip. The distorted
# particle balance takes its square root, so (1 - alpha)^0.5 is already present before any
# hindrance is applied. The residual that recovers Ishii's (1 - alpha)^0.75 relative velocity is
# therefore 0.25, the same value hibiki-annulus.i carries and for the same reason.
COMMON="jf=0.491 dp=0.0025 drag_model=distorted-particle swarm_exponent=0.25"
STEPS="Executioner/num_steps=2000 Executioner/dt=0.02 Outputs/checkpoint=false"
LTOL="Executioner/momentum_l_tol=1e-10 Executioner/pressure_l_tol=1e-10 Executioner/active_scalar_l_tol=1e-10"

run () {  # run <model> <alpha> <jg> <ai_in>
  # The density time derivative is selected with the density: non-zero only when rho_d follows
  # the solved pressure.
  if [ "$1" = "solved" ]; then DRHO="p_dot drho_g_dt_solved"; else DRHO="drho_g_dt_zero"; fi
  $EXE -i $IN "FunctorMaterials/active=rho_g_$1 $DRHO c_d j_axial alpha_j alpha_v_gj" \
       alpha_in=$2 jg=$3 ai_in=$4 $COMMON $STEPS $LTOL \
       Outputs/file_base=v5_$1_$2 > v5_$1_$2.log 2>&1
}

# One density model at a time, its four void fractions in parallel. Launching all twelve at once
# is slower in aggregate than four: the runs are memory bandwidth bound rather than core bound, and
# oversubscribing drops the total throughput below what a smaller batch sustains.
for m in solved imposed constant; do
  run $m 0.05 0.0275 93.2  &
  run $m 0.10 0.0556 176.3 &
  run $m 0.20 0.129  321.0 &
  run $m 0.25 0.170  406.3 &
  wait
  echo "batch $m done"
done
echo "SWEEP DONE"
