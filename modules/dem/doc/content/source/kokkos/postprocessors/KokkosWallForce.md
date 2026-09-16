# KokkosWallForce

!syntax description /Postprocessors/KokkosWallForce

Reports a component, or the magnitude, of the total contact force the particles of a
[KokkosParticleCloud.md] exert on one of its walls: an analytic wall by its index in
`wall_points` (`wall`), or a sideset wall by boundary name (`boundary`). The force is the sum over
the particles' wall contacts of the normal and tangential forces from the current contact
histories, negated, at the end of the step, reduced across ranks. The `stack_force` test reads
the weight of a settled stack on its floor.

!syntax parameters /Postprocessors/KokkosWallForce

!syntax inputs /Postprocessors/KokkosWallForce
