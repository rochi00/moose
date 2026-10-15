# KokkosLinearFVFunctorDirichletBC

!if! function=hasCapability('kokkos')

This object provides Dirichlet boundary value data for Kokkos linear finite volume kernels. The
boundary value is supplied through a Kokkos-compatible `functor`.

For kernels that need boundary normal gradient data, this object supplies the affine relation

!equation
\nabla_n u = -\frac{1}{d_{Cf}} u_C + \frac{g}{d_{Cf}}

where $u_C$ is the cell value, $g$ is the prescribed boundary value, and $d_{Cf}$ is the distance
from the cell center to the face along the face normal, $|\mathbf{d}_{Cf} \cdot \mathbf{n}|$, as in
[LinearFVAdvectionDiffusionFunctorDirichletBC.md]. Kernels apply their own physics coefficients when consuming
this relation.

The `functor` is a number, a [KokkosConstantFunction.md], a [KokkosParsedFunction.md], a linear
finite volume variable or a [KokkosParsedFunctorMaterial.md] property; other function types wait
on [relocatable device code (RDC)](KokkosFunctions/index.md#kokkos_rdc) support. As the host
object, it evaluates the functor without a side, so a function at the face centroid also on block
interfaces, and a functor material there as the interpolation of its cell values.

## Example Syntax

!listing test/tests/kokkos/linearfvkernels/diffusion/kokkos_diffusion-2d.i start=[dir] end=[] include-end=true

!syntax parameters /LinearFVBCs/KokkosLinearFVFunctorDirichletBC

!syntax inputs /LinearFVBCs/KokkosLinearFVFunctorDirichletBC

!syntax children /LinearFVBCs/KokkosLinearFVFunctorDirichletBC

!if-end!

!else
!include kokkos/kokkos_warning.md
