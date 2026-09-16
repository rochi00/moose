# KokkosParticleProfile

!syntax description /VectorPostprocessors/KokkosParticleProfile

Coarse-grains the particles of a [KokkosParticleCloud.md] into `num_bins` slabs of equal width
along a `direction`, over the mesh bounding box or the given `lower` to `upper`, each particle
counted in the slab of its center. Per slab it reports the number of particles, the solid
fraction (the particles' volume $\frac{4}{3}\pi r^3$ over the slab's volume, its width times the
bounding box's extent in the other directions of the mesh dimension), the mean velocity, the
granular temperature $T = m \langle |v - \bar v|^2 \rangle / 3$, and the stress tensor

!equation
\sigma = -\frac{1}{V} \left( \sum_i m_i v'_i \otimes v'_i + \sum_\text{pairs} r_{ij} \otimes f_{ij}
+ \sum_\text{walls} r_{iw} \otimes f_{iw} \right)

with the fluctuation velocities $v'$ about the slab mean, half of each pair's virial going to the
slab of each particle, and the wall term taken from the contact point to the center. The sums are
reduced across ranks. The `profile` test checks the kinetic part on the initial gas of the
`periodic_gas` case and `stack_profile` the contact part on the settled Hertz stack.

!syntax parameters /VectorPostprocessors/KokkosParticleProfile

!syntax inputs /VectorPostprocessors/KokkosParticleProfile
