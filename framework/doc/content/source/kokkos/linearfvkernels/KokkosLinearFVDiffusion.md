# KokkosLinearFVDiffusion

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVDiffusion.md]. As there, the nonorthogonal correction
(`use_nonorthogonal_correction`, on by default) splits the normal gradient on a face into an
implicit difference of the cell values along the cell-to-cell vector and an explicit correction
from the face gradient, linearly interpolated from the cell gradients of the variable's gradient
method. Boundary conditions that compute the boundary value from the cell value, such as
[KokkosLinearFVFunctorDirichletBC.md], get the correction on their faces as well, from the cell
gradient. The correction needs the cell gradients, which Green-Gauss gradient methods do not
support in spherical coordinates; on orthogonal meshes, where the correction vanishes, disable it.

## Example Syntax

!listing test/tests/kokkos/linearfvkernels/diffusion/kokkos_diffusion-2d.i start=[diffusion] end=[] include-end=true

!syntax parameters /LinearFVKernels/KokkosLinearFVDiffusion

!syntax inputs /LinearFVKernels/KokkosLinearFVDiffusion

!syntax children /LinearFVKernels/KokkosLinearFVDiffusion

!if-end!

!else
!include kokkos/kokkos_warning.md
