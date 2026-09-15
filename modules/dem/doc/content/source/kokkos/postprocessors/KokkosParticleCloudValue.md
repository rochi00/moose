# KokkosParticleCloudValue

!syntax description /Postprocessors/KokkosParticleCloudValue

Reports one of the global scalars of a [KokkosParticleCloud.md]: the number of particles, the
number lost, the number migrated between ranks in the last step, the maximum number of face hops
any particle took in the last step, the translational and rotational kinetic energies, or the
magnitude of the total spin angular momentum. All values are reduced across ranks.

!syntax parameters /Postprocessors/KokkosParticleCloudValue

!syntax inputs /Postprocessors/KokkosParticleCloudValue

!syntax children /Postprocessors/KokkosParticleCloudValue
