# KokkosFVGreenGaussGradient

!if! function=hasCapability('kokkos')

This is the Kokkos version of [FVGreenGaussGradient.md]. It writes the same gradients into the
same storage, so both host and Kokkos linear finite volume objects can read them. On boundary
faces of a variable the face value comes from the value relation of the Kokkos linear finite volume
boundary condition on that face, or from the cell value where there is none. Each cell sums its
faces in a fixed order, so the gradients are identical in every run. The `venkatakrishnan`
limiter is applied on the device as well, with the same formula as the host limiter.

[FVGreenGaussGradient.md] uses the same Kokkos computation, limiter included, for variables whose
linear finite volume kernels and boundary conditions are all Kokkos objects.

## Example Syntax

!listing test/tests/kokkos/fvgradientmethods/green-gauss.i block=FVGradientMethods

!syntax parameters /FVGradientMethods/KokkosFVGreenGaussGradient

!syntax inputs /FVGradientMethods/KokkosFVGreenGaussGradient

!syntax children /FVGradientMethods/KokkosFVGreenGaussGradient

!if-end!

!else
!include kokkos/kokkos_warning.md
