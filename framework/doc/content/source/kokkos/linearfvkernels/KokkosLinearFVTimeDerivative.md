# KokkosLinearFVTimeDerivative

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVTimeDerivative.md] for a constant multiplier
([!param](/LinearFVKernels/KokkosLinearFVTimeDerivative/factor)). It supports the implicit Euler
and BDF2 time integrators.

## Example Syntax

!listing test/tests/kokkos/linearfvkernels/time-derivative/time-derivative.i start=[time] end=[] include-end=true

!syntax parameters /LinearFVKernels/KokkosLinearFVTimeDerivative

!syntax inputs /LinearFVKernels/KokkosLinearFVTimeDerivative

!syntax children /LinearFVKernels/KokkosLinearFVTimeDerivative

!if-end!

!else
!include kokkos/kokkos_warning.md
