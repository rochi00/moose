# KokkosParticleCloudValue

!syntax description /Postprocessors/KokkosParticleCloudValue

Reports one of the global scalars of a [KokkosParticleCloud.md]: the number of particles, the
number that exited the mesh (cumulative), the number left unresolved by the face walk in the last
step, the number migrated between ranks in the last step, the maximum number of face hops any
particle took in the last step, the number of neighbor pairs (at the last build), neighbor-list
builds, ghost particles, and ghost forwards, the number of overlapping pairs and the coordination
number (twice that over the number of particles), the number of contact histories held (summed
over ranks), the load imbalance (the largest local particle count over the mean), the translational, rotational, spring (pairs and walls), and total energies,
the magnitudes of the total linear and spin angular momenta, or the virial pressure

!equation
P = \frac{1}{d V} \left( 2 K + \sum_{i<j} r_{ij} \cdot f_{ij} \right)

over the overlapping pairs, with $V$ the volume of the mesh bounding box in the mesh dimension
$d$ and $K$ the translational kinetic energy; or the extent of the bed, the largest and smallest
coordinate of a particle surface in each direction (`max_x` .. `min_z`), and the center of mass
(`center_of_mass_x` ..). All values are reduced across ranks.

!syntax parameters /Postprocessors/KokkosParticleCloudValue

!syntax inputs /Postprocessors/KokkosParticleCloudValue

!syntax children /Postprocessors/KokkosParticleCloudValue
