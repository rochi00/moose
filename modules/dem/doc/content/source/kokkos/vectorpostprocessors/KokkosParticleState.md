# KokkosParticleState

!syntax description /VectorPostprocessors/KokkosParticleState

Outputs the global ID, position, velocity, and body-frame angular velocity of every particle of a
[KokkosParticleCloud.md], gathered to the root rank and sorted by global ID so the output does not
depend on the partitioning.

!syntax parameters /VectorPostprocessors/KokkosParticleState

!syntax inputs /VectorPostprocessors/KokkosParticleState

!syntax children /VectorPostprocessors/KokkosParticleState
