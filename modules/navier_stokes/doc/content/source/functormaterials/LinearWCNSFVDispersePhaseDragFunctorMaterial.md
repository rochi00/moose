# LinearWCNSFVDispersePhaseDragFunctorMaterial

Computes the linear drag function between the dispersed and the continuous phase of a two-phase
mixture, for the linear finite volume discretization. It is computed in plain reals rather than
automatic differentiation types, since the linear finite volume discretization assembles a matrix
directly.

The drag function is that of Schiller and Naumann below the transition, and Newton's regime above
it:

\begin{equation}
f_{drag} =
\begin{cases}
1 + 0.15 Re_p^{0.687}, & Re_p \leq 1000 \\
0.0183 Re_p,           & Re_p > 1000
\end{cases}
\end{equation}

Both a scalar property and a vector property, suffixed `_vec`, are declared. The vector form is what
the momentum friction kernels consume.

with the particle Reynolds number formed from the slip velocity and the continuous phase
properties, which is its definition, see [!cite](manninen1996mixture) equation (39):

!equation
Re_p = \frac{\rho_c d_d \|\bm{u}_{slip,d}\|}{\mu_c}

Accordingly the `u`, `v` and `w` parameters take the *slip* velocity components, and `rho` and `mu`
take the *continuous phase* density and viscosity, not the mixture ones.

!alert note title=Why this does not close a loop
The slip velocity closure also needs a drag function, so taking the slip velocity as an input here
looks circular. It is not, because
[LinearWCNSFV2PSlipVelocityFunctorMaterial.md] solves the drag correlation together with its own
force balance rather than consuming a drag functor. This object is therefore downstream of the slip
velocity and nothing depends on it in turn, except the momentum friction kernels that consume the
`_vec` form. Combining `use_dispersed_phase_drag_model` with a prescribed
`slip_linear_friction_name` on the Physics is rejected for the same reason.

!syntax parameters /FunctorMaterials/LinearWCNSFVDispersePhaseDragFunctorMaterial

!syntax inputs /FunctorMaterials/LinearWCNSFVDispersePhaseDragFunctorMaterial

!syntax children /FunctorMaterials/LinearWCNSFVDispersePhaseDragFunctorMaterial
