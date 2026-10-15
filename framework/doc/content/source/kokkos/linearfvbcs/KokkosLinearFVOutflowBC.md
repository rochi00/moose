# KokkosLinearFVOutflowBC

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVAdvectionDiffusionOutflowBC.md]. The boundary value is extrapolated from the cell value, with the cell gradient along the cell-to-face vector if `use_two_term_expansion` is set, and there is no diffusive flux.

## Example Syntax

!listing test/tests/kokkos/linearfvbcs/framework-bcs.i start=[outflow] end=[] include-end=true

!syntax parameters /LinearFVBCs/KokkosLinearFVOutflowBC

!syntax inputs /LinearFVBCs/KokkosLinearFVOutflowBC

!syntax children /LinearFVBCs/KokkosLinearFVOutflowBC

!if-end!

!else
!include kokkos/kokkos_warning.md
