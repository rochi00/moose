# KokkosLinearFVFunctorRobinBC

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVAdvectionDiffusionFunctorRobinBC.md]. It imposes $\alpha \nabla u \cdot \mathbf{n} + \beta u = \gamma$ with the same boundary value and normal gradient relations as the host object, including the cell gradient along the tangential part of the cell-to-face vector, and the projections on the normal of the face information. The coefficients `alpha`, `beta` and `gamma` are numbers, Kokkos functions, linear finite volume variables or [KokkosParsedFunctorMaterial.md] properties, evaluated from the side of the variable.

## Example Syntax

!listing test/tests/kokkos/linearfvbcs/framework-bcs.i start=[robin] end=[] include-end=true

!syntax parameters /LinearFVBCs/KokkosLinearFVFunctorRobinBC

!syntax inputs /LinearFVBCs/KokkosLinearFVFunctorRobinBC

!syntax children /LinearFVBCs/KokkosLinearFVFunctorRobinBC

!if-end!

!else
!include kokkos/kokkos_warning.md
