# KokkosParticleCloudValue

!syntax description /Postprocessors/KokkosParticleCloudValue

Reports one of the global scalars of a [KokkosParticleCloud.md]: the number of particles, the
number that exited the mesh (cumulative), the number left unresolved by the face walk in the last
step, the number migrated between ranks in the last step, the maximum number of face hops any
particle took in the last step, the number of neighbor pairs (at the last build), neighbor-list
builds, ghost particles, and ghost forwards, the translational and rotational kinetic energies, or
the magnitude of the total spin angular momentum. All values are reduced across ranks.

!syntax parameters /Postprocessors/KokkosParticleCloudValue

!syntax inputs /Postprocessors/KokkosParticleCloudValue

!syntax children /Postprocessors/KokkosParticleCloudValue
