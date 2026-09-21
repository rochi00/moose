# Drift-flux model: linear FV segregated vs nonlinear FV Newton

Two cases run with both discretizations of the mixture drift-flux model, to compare the answers
and the cost of getting them. The inputs share their model definition through `!include`, so
each pair differs only in the solver blocks.

## tjunction

Three-dimensional T-junction at Re = 100 with a light dispersed phase and gravity across the
main channel, 576 cells at the base resolution (`Mesh/uniform_refine=1` for 4608).

- `simple.i`: linear FV, SIMPLE steady solve (no gravity effect on the phase distribution is
  visible without a transient, see `pimple.i`)
- `pimple.i`: linear FV, PIMPLE time march to steady state
- `newton.i`: nonlinear FV, Newton time march to steady state, same time step

## pipe

Vertical bubbly air-water flow in the 50.8 mm pipe of Hibiki, Ishii and Xiao (2001),
axisymmetric, 1920 cells, 400 steps of 0.02 s. Postprocessors form the one-dimensional
drift-flux parameters C0 and Vgj at the three measuring stations of the experiment.

- `hibiki-50mm-pipe.i`: linear FV, PIMPLE, algebraic multigrid on every system
- `newton-pipe.i`: nonlinear FV, Newton, direct factorization

## Running

    ./run.sh tjunction/pimple.i
    NP=4 ./run.sh pipe/hibiki-50mm-pipe.i Mesh/uniform_refine=1
    NP=4 ./run.sh pipe/newton-pipe.i MUMPS Mesh/uniform_refine=1

`run.sh` appends one line per run to `results.log` with the wall time, number of steps, number
of failed solves and, for the pipe, C0 and Vgj at z/D = 54. `MUMPS` on the command line hands
the Newton factorization to MUMPS, which is needed for a direct solve on more than one rank.
