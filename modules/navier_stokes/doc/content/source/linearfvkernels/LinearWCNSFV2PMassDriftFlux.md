# LinearWCNSFV2PMassDriftFlux

Adds the dilatation produced by the relative motion of the phases to the pressure equation of a two
phase mixture.

## The inconsistency this removes

The mixture momentum equation is written for the mass averaged velocity $\bm{u}_m$, whose continuity
equation is

\begin{equation}
\frac{\partial \rho_m}{\partial t} + \nabla \cdot \left(\rho_m \bm{u}_m\right) = 0 .
\end{equation}

$\bm{u}_m$ is therefore *not* solenoidal wherever the mixture density varies, and the mixture
density varies with the phase fraction. Constraining $\nabla\cdot\bm{u}_m = 0$, as the single phase
pressure correction of [SIMPLE.md] does, imposes the constraint belonging to the volume averaged
velocity on a field that is mass averaged, and omits the dilatation an evolving phase fraction
produces.

## The form assembled

With constant phase densities and no phase change the two phase volume equations sum to

\begin{equation}
\nabla \cdot \bm{j} = 0 , \qquad \bm{j} = \alpha \bm{u}_d + \left(1-\alpha\right)\bm{u}_c ,
\end{equation}

which is an exact restatement of the mixture continuity equation above rather than an additional
assumption. Using the identity $\bm{j} = \bm{u}_m + \left(\alpha - c_d\right)\bm{u}_s$ with
$c_d = \alpha\rho_d/\rho_m$ the mass fraction of the dispersed phase, the constraint the pressure
equation should impose is

\begin{equation}
\nabla \cdot \bm{u}_m = -\nabla \cdot \left[\left(\alpha - c_d\right)\bm{u}_s\right] .
\end{equation}

This object assembles the right hand side, face by face and with the same sign convention as
[LinearFVDivergence.md] uses for the divergence of the predicted flux.

Because it is an algebraic identity it carries no time derivative. It is therefore equally valid in
a steady solve, where $\bm{u}_m$ is still not solenoidal because $\rho_m$ varies in space, and it
raises no question of consistency between the discrete time derivative here and the one in the phase
equation.

The contribution is withheld on impermeable boundaries, where the phases cannot separate and the
relative motion produces no dilatation.

## Verification

The factor $\alpha - c_d$ vanishes identically when the two phase densities are equal, since
$c_d = \alpha$ then and the volume averaged and mass averaged velocities coincide. The
`mass_drift_flux_vanishes` test asserts exactly that: run with equal phase densities and the term
active, it reproduces to the last digit a gold file generated with the term switched off. This is an
analytic property rather than a tolerance.

!alert note
The magnitude of the term is set by the *divergence* of $\left(\alpha - c_d\right)\bm{u}_s$, not by
its magnitude. In a developed channel or pipe that field is nearly uniform along the flow, so the
correction is small there even when the slip is not; it matters where the phase fraction varies
sharply.

!syntax parameters /LinearFVKernels/LinearWCNSFV2PMassDriftFlux

!syntax inputs /LinearFVKernels/LinearWCNSFV2PMassDriftFlux

!syntax children /LinearFVKernels/LinearWCNSFV2PMassDriftFlux
