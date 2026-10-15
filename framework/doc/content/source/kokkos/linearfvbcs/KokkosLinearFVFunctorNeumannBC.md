# KokkosLinearFVFunctorNeumannBC

!if! function=hasCapability('kokkos')

This is the Kokkos version of [LinearFVAdvectionDiffusionFunctorNeumannBC.md]. It prescribes the
outward diffusive flux $q$, supplied as a number, a [KokkosConstantFunction.md], a
[KokkosParsedFunction.md], a linear finite volume variable or a [KokkosParsedFunctorMaterial.md]
property via the [!param](/LinearFVBCs/KokkosLinearFVFunctorNeumannBC/functor) parameter. The
flux and the diffusion coefficient are evaluated from the variable's side of the face, as on the
host. Diffusion kernels add
$q$ times the face area to the right-hand side, without their own diffusion coefficient.

For kernels that need a boundary value, such as advection, and for the Green-Gauss gradient, it
extrapolates the cell value $u_C$:

!equation
u_b = u_C + \frac{q}{k} \left|\mathbf{d}_{Cf} \cdot \mathbf{n}\right| + \nabla u_C \cdot
\left(\mathbf{d}_{Cf} - (\mathbf{d}_{Cf} \cdot \mathbf{n}) \mathbf{n}\right)

where $k$ is the diffusion coefficient given by
[!param](/LinearFVBCs/KokkosLinearFVFunctorNeumannBC/diffusion_coeff), $\mathbf{d}_{Cf}$ the
vector from the cell centroid to the face centroid, and $\mathbf{n}$ the outward unit normal. The
last term corrects for non-orthogonal meshes using the cell gradient of the variable, so the
variable's [!param](/Variables/MooseLinearVariableFVReal/gradient_method) should be a Kokkos
method such as [KokkosFVGreenGaussGradient.md], which uses the Kokkos boundary conditions.

## Example Syntax

!listing test/tests/kokkos/linearfvkernels/diffusion/kokkos_diffusion-2d-neumann.i block=LinearFVBCs

!syntax parameters /LinearFVBCs/KokkosLinearFVFunctorNeumannBC

!syntax inputs /LinearFVBCs/KokkosLinearFVFunctorNeumannBC

!syntax children /LinearFVBCs/KokkosLinearFVFunctorNeumannBC

!if-end!

!else
!include kokkos/kokkos_warning.md
