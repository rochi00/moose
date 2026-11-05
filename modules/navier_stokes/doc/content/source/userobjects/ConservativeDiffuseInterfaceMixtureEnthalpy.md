# ConservativeDiffuseInterfaceMixtureEnthalpy

!syntax description /UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy

## Description

`ConservativeDiffuseInterfaceMixtureEnthalpy` relates the mixture enthalpy density $\rho h$ of a
two-phase diffuse-interface flow to its temperature $T$:

\begin{equation}
\rho h = \alpha \left\{ f_l \rho_l \left[ h_l(T) + L \right] + (1 - f_l) \rho_s h_s(T) \right\}
       + (1 - \alpha) \rho_g h_g(T),
\end{equation}

where $\alpha$ is the liquid
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/volume_fraction), $\rho_l$ and
$\rho_g$ are the constant
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/liquid_density) and
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/gas_density), and $h_l$ and
$h_g$ are the specific enthalpies of
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/liquid_fluid_properties) and
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/gas_fluid_properties) at the
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/thermodynamic_pressure).
The liquid phase melts and solidifies when the
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/latent_heat) $L$ is positive;
$f_l$ is the volume fraction of its liquid, and $\rho_s$ and $h_s$ the density and enthalpy of
its solid (see [#phase-change]). Without latent heat the braces hold $\rho_l h_l(T)$.

It does three things:

- At the start of the simulation, it sets
  [!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/enthalpy_variable) from the
  initial [!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/temperature) and
  volume fraction. A recovered simulation restores $\rho h$ from the checkpoint instead.
- After each energy solve, it recovers the temperature in every cell by Newton iteration on the
  relation above, starting from the previous temperature. It stops with an error if a cell does
  not converge.
- It provides the energy kernels with the linearization of $\rho h$ about the lagged state of a
  cell, $T = T^* + \frac{dT}{d(\rho h)} \left( \rho h - (\rho h)^* \right)$, and the matching
  derivatives $dh_k / d(\rho h)$ of the phase enthalpies. Without phase change,
  $dT/d(\rho h) = 1/C$ with the mixture heat capacity $C = \sum_k \alpha_k \rho_k c_{p,k}(T^*)$.
  It also provides the cell thermal conductivity, the series mix $1/k = \sum_k \alpha_k / k_k$.

The series mix makes the cell resistances across a diffuse interface add up to those of the two
sharp layers, however many cells the interface spreads over.

## Phase change id=phase-change

With latent heat, the liquid phase is fully solid below the
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/solidus_temperature) $T_s$ and
fully liquid above the
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/liquidus_temperature) $T_l$.
Between them its liquid fraction $f_l(T)$ is linear or the smoothstep $3 s^2 - 2 s^3$ of
$s = (T - T_s) / (T_l - T_s)$, as chosen by
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/liquid_fraction_shape). Equal
solidus and liquidus temperatures describe melting at a single temperature $T_m$, where a cell
holds any liquid fraction between 0 and 1. The liquid fraction is therefore stored in the
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/liquid_fraction) variable and
recovered from $\rho h$ together with the temperature, rather than computed from $T$.

The recovery compares $\rho h$ with the cell's fully solid enthalpy $H_s$ at $T_s$ and fully liquid
enthalpy $H_l$ at $T_l$:

- below $H_s$ the cell is solid, $f_l = 0$ and $T \le T_s$;
- above $H_l$ it is liquid, $f_l = 1$ and $T \ge T_l$;
- in between it melts: over a finite interval $T$ solves the relation above with $f_l = f_l(T)$;
  at a single melting temperature $T = T_m$ and $f_l = (\rho h - H_s) / (H_l - H_s)$.

The enthalpy is smooth in $T$ on each of these branches. The Newton iteration keeps a bracket of
the root and bisects when a step leaves it, which is needed where the slope of the melting branch
is far from its secant. In the linearization, $dT/d(\rho h) = 1/(C + \alpha \rho_l L f_l'(T^*))$
over a finite interval and, for a cell at $T_m$, that of a latent heat spread over 0.01 K, so that
the cell stays at $T_m$ while it changes phase. The width keeps the derivative finite: held exactly
at $T_m$, a cell whose latent content is small against the heat it exchanges in a step, such as a
gas cell with a trace of liquid, would give up far more enthalpy than it holds. The latent heat therefore cannot be skipped by a time step that carries a cell across the
melting interval, and $\rho h$ stays conserved.

A cell takes part in the phase change only when its latent heat $\alpha \rho_l L$ is resolvable
against the round-off of $\rho h$, which carries the absolute phase enthalpies: it must exceed
$\sqrt{\epsilon} \, |\rho h|$, with $\epsilon$ the machine precision. Cells with less liquid, in
particular all gas cells, keep their liquid fraction unchanged.

A boundary or inflow temperature implies the liquid fraction $f_l(T)$, with a single melting
temperature counted as liquid.

The solid takes its heat capacity and thermal conductivity from
[!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/solid_fluid_properties) and its
density from [!param](/UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy/solid_density), or
shares the liquid's when these are omitted. The solid enthalpy is shifted to equal the liquid's at
$T_s$, so that $L$ is the whole enthalpy difference between the phases there whatever reference
states the two fluid properties objects use. Within a cell the liquid and solid conduct in
parallel, $k = f_l k_l + (1 - f_l) k_s$, before the series mix with the gas; in a mushy cell they
are the same material side by side. The latent capacity of the linearization and of the Newton
iteration over a finite interval is $\alpha \, d(\rho h)/df_l = \alpha [\rho_l (h_l + L) - \rho_s h_s]$,
which is $\alpha \rho_l L$ at $T_s$ for equal densities and grows away from it with a separate
solid heat capacity.

The liquid mass flux of the energy equation advects the enthalpy of the liquid, $h_l + L$, since
the solid is held at rest. A solid density different from the liquid's changes the volume of the
liquid phase as it changes phase: a freezing cell expels liquid, which carries the latent heat away
with it. This object only accounts for the enthalpy each unit of volume then holds; the matching
volume source of the flow, $\nabla \cdot \mathbf{u} = \alpha \, (1 - \rho_s / \rho_l) \,
\partial f_s / \partial t$ with $f_s = 1 - f_l$ [!cite](lyu2021hybrid), is added to the pressure
and volume-fraction equations in the input.

This object is created by
[WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics](WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics.md),
which also creates the kernels and boundary conditions that use it.

## Example Input Syntax

In this example, the heat transfer physics creates this object for the conserved enthalpy
`rho_h` and the temperature `T` of a liquid layer under a gas layer. The two phases have different
conductivities and heat capacities that vary with temperature, given by the fluid properties
objects `fp_l` and `fp_g`.

!listing modules/navier_stokes/test/tests/finite_volume/two_phase/diffuse_interface/two_layer_conduction.i block=Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated

!syntax parameters /UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy

!syntax inputs /UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy

!syntax children /UserObjects/ConservativeDiffuseInterfaceMixtureEnthalpy
