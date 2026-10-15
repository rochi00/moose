# KokkosLinearFVExtrapolatedBC

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVAdvectionDiffusionExtrapolatedBC.md]. The boundary value is extrapolated from the cell value, with the cell gradient along the cell-to-face vector if `use_two_term_expansion` is set, and the normal gradient is the difference of the boundary value, computed with the current cell value, and the cell value over the normal distance to the face.

## Example Syntax

!listing test/tests/kokkos/linearfvbcs/framework-bcs.i start=[extrapolated] end=[] include-end=true

!syntax parameters /LinearFVBCs/KokkosLinearFVExtrapolatedBC

!syntax inputs /LinearFVBCs/KokkosLinearFVExtrapolatedBC

!syntax children /LinearFVBCs/KokkosLinearFVExtrapolatedBC

!if-end!

!else
!include kokkos/kokkos_warning.md
