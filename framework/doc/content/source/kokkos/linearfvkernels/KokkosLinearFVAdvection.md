# KokkosLinearFVAdvection

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVAdvection.md]. The advected face values on internal faces
come from the `FVInterpolationMethods` object named by `advected_interp_method_name`, which must be
an [FVAdvectedUpwind.md], [FVGeometricAverage.md], [FVAdvectedMinmodWeightBased.md],
[FVAdvectedVanLeerWeightBased.md] or [FVAdvectedMUSCLDeferredCorrection.md]; the kernel reads its
parameters and evaluates it on the device. The methods that need cell gradients read them from the
gradient method the interpolation method names.

!alert note
`KokkosLinearFVAdvection` supports constant velocities only.

## Example Syntax

!listing test/tests/kokkos/linearfvkernels/advection-interpolation/advection-interpolation.i block=FVInterpolationMethods

!listing test/tests/kokkos/linearfvkernels/advection-interpolation/advection-interpolation.i start=[advection] end=[] include-end=true

!syntax parameters /LinearFVKernels/KokkosLinearFVAdvection

!syntax inputs /LinearFVKernels/KokkosLinearFVAdvection

!syntax children /LinearFVKernels/KokkosLinearFVAdvection

!if-end!

!else
!include kokkos/kokkos_warning.md
