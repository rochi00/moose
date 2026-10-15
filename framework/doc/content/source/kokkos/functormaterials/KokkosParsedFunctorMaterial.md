# KokkosParsedFunctorMaterial

!if! function=hasCapability('kokkos')

This is the Kokkos version of [ParsedFunctorMaterial.md], with the same parameters. It defines a
property of the Kokkos linear finite volume kernels and boundary conditions, such as a diffusion
coefficient, from a parsed expression of linear finite volume variables (solver or auxiliary),
numbers, space (`x`, `y`, `z`) and time (`t`); `pi` and `e` are defined as in
[ParsedFunctorMaterial.md]. The objects refer to it by its `property_name`, and evaluate it on the
device where the host objects evaluate a functor material:

- on a cell, the expression of the cell values of the variables at the cell centroid;
- on a face with the object's variable on both sides, and on a face between blocks evaluated
  without a side, the linear interpolation of the cell values of the property;
- on a face evaluated from the variable's side, such as a boundary face for a diffusion
  coefficient, the expression of the face values of the variables at the face centroid. The face
  value of a variable is the linear interpolation of its cell values where it lives on both sides,
  the boundary value of its Kokkos boundary condition, or its cell value.

The values are recomputed before every assembly that uses them, as they depend on the solution.
Unlike [ParsedFunctorMaterial.md], it defines no host functor, so host objects cannot use the
property, and the functors in the expression are limited to linear finite volume variables and
numbers; the dependence on space and time goes in the expression itself.

## Example Syntax

!listing test/tests/kokkos/linearfvkernels/properties/variable-properties.i block=FunctorMaterials

!listing test/tests/kokkos/linearfvkernels/properties/variable-properties.i start=[diffusion_u] end=[] include-end=true

!syntax parameters /FunctorMaterials/KokkosParsedFunctorMaterial

!syntax inputs /FunctorMaterials/KokkosParsedFunctorMaterial

!syntax children /FunctorMaterials/KokkosParsedFunctorMaterial

!if-end!

!else
!include kokkos/kokkos_warning.md
