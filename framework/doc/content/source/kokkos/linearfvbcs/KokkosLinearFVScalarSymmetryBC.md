# KokkosLinearFVScalarSymmetryBC

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVAdvectionDiffusionScalarSymmetryBC.md]. There is no normal gradient, and the boundary value is the cell value, with the cell gradient along the tangential part of the cell-to-face vector if `use_two_term_expansion` is set.

## Example Syntax

!listing test/tests/kokkos/linearfvbcs/framework-bcs.i start=[symmetry] end=[] include-end=true

!syntax parameters /LinearFVBCs/KokkosLinearFVScalarSymmetryBC

!syntax inputs /LinearFVBCs/KokkosLinearFVScalarSymmetryBC

!syntax children /LinearFVBCs/KokkosLinearFVScalarSymmetryBC

!if-end!

!else
!include kokkos/kokkos_warning.md
