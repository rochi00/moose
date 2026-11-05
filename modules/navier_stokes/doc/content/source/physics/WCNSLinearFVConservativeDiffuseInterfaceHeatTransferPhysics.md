# WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics

!syntax description /Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics

## Equation

`WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics` adds the energy equation to a
two-phase diffuse-interface flow solved with the linear finite volume reduced-pressure path
described in the [diffuse-interface dam break algorithm](diffuse_interface_dam_break_algorithm.md)
note. The solved unknown is the mixture enthalpy
density $\rho h$, which is conserved across the interface:

\begin{equation}
\frac{\partial (\rho h)}{\partial t}
+ \nabla \cdot \left( \sum_k \dot{m}_k h_k(T) \right)
- \nabla \cdot \left( k \nabla T \right) = 0,
\qquad
\rho h = \sum_k \alpha_k \rho_k h_k(T) + \alpha \rho_l L f_l,
\end{equation}

where:

- $k \in \{l, g\}$ are the liquid and gas phases, with volume fractions $\alpha_l = \alpha$ and
  $\alpha_g = 1 - \alpha$
- $\rho_k$ are the constant phase densities of the coupled volume-fraction physics
- $h_k(T)$ are the phase specific enthalpies at the
  [!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/thermodynamic_pressure),
  which may depend nonlinearly on the temperature
- $\dot{m}_k$ are the phase mass fluxes through a face. The liquid flux is
  $\dot{m}_l = \rho_l F_\alpha$ and the gas flux is $\dot{m}_g = F_\rho - \dot{m}_l$, where
  $F_\alpha$ and $F_\rho$ are the limited volume-fraction flux `alpha_phi_limited` and the mass
  flux `rho_phi` that advanced the volume fraction, so energy is carried with exactly the mass
  each phase moved
- $k$ is the mixture thermal conductivity, $1/k = \sum_k \alpha_k / k_k$, so the phases conduct
  in series across the interface
- $L$ is the latent heat and $f_l$ the liquid fraction of the liquid phase, nonzero only when it
  melts and solidifies (see [#phase-change]); the liquid mass flux then carries $h_l + L$, the
  enthalpy of the liquid, as the solid is at rest

The temperature $T$ is an auxiliary variable recovered from $\rho h$ in each cell after every
energy solve, by Newton iteration on the definition of $\rho h$ above.

The objects created are:

- [LinearFVTimeDerivative](LinearFVTimeDerivative.md) for the time derivative of $\rho h$
- [LinearFVConservativeDiffuseInterfaceEnthalpyAdvection](LinearFVConservativeDiffuseInterfaceEnthalpyAdvection.md)
  for the phase-wise advection
- [LinearFVConservativeDiffuseInterfaceEnthalpyDiffusion](LinearFVConservativeDiffuseInterfaceEnthalpyDiffusion.md)
  for heat conduction
- [LinearFVConservativeDiffuseInterfaceEnthalpyBC](LinearFVConservativeDiffuseInterfaceEnthalpyBC.md)
  on every inlet, wall and outlet of the coupled flow
- [ConservativeDiffuseInterfaceMixtureEnthalpy](ConservativeDiffuseInterfaceMixtureEnthalpy.md),
  which relates $\rho h$ to $T$, sets the initial $\rho h$ and recovers $T$ (and $f_l$)

## Phase properties

Each phase takes its specific enthalpy, heat capacity and thermal conductivity from a
[SinglePhaseFluidProperties](SinglePhaseFluidProperties.md) object, given by
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/liquid_fluid_properties)
and
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/gas_fluid_properties).
They are evaluated at the constant
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/thermodynamic_pressure).
The densities are not taken from the fluid properties: they are the constant phase densities of
the coupled volume-fraction physics. A
[TemperaturePressureFunctionFluidProperties](TemperaturePressureFunctionFluidProperties.md) object
gives a heat capacity that varies with temperature.

## Melting and solidification id=phase-change

A positive
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/latent_heat)
lets the liquid phase melt and solidify between the
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/solidus_temperature)
and the
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/liquidus_temperature),
or at a single melting temperature when the two are equal. The volume fraction $\alpha$ then
counts liquid and solid together, so its transport is unchanged, and the solid fraction of a cell
is $\alpha (1 - f_l)$. The physics creates the liquid fraction as the auxiliary variable
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/liquid_fraction_variable),
which [ConservativeDiffuseInterfaceMixtureEnthalpy](ConservativeDiffuseInterfaceMixtureEnthalpy.md)
recovers from $\rho h$ with the temperature. The solid takes its heat capacity and thermal
conductivity from
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/solid_fluid_properties)
and its density from
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/solid_density),
or the liquid's when these are omitted. A different solid density changes the volume of the
liquid phase as it freezes or melts; the volume source that goes with it is added to the flow in
the input, as described in
[ConservativeDiffuseInterfaceMixtureEnthalpy](ConservativeDiffuseInterfaceMixtureEnthalpy.md).

The energy equation does not hold the solid still. A momentum sink that does, such as a
Carman-Kozeny drag on the effective liquid fraction $1 - \alpha (1 - f_l)$, is added in the input,
for example with
[INSFVMushyPorousFrictionFunctorMaterial](INSFVMushyPorousFrictionFunctorMaterial.md) and a
[LinearFVReaction](LinearFVReaction.md) on each velocity component.

## Boussinesq buoyancy

The phase densities in the continuity, momentum and energy equations stay constant. A
temperature-dependent density can enter the gravity force only, through the Boussinesq density

\begin{equation}
\rho_b = \sum_k \alpha_k \rho_k \left( 1 - \beta_k \left( T - T_{ref} \right) \right),
\end{equation}

where $\beta_k$ are
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/liquid_thermal_expansion)
and
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/gas_thermal_expansion),
and $T_{ref}$ is
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/reference_temperature),
which is required when either expansion coefficient is nonzero. The density is defined as a
[ParsedFunctorMaterial](ParsedFunctorMaterial.md) named by
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/buoyancy_density_name).
It takes effect only when that name is also passed as the `buoyancy_density` of the coupled
`ConservativeDiffuseInterfaceFlowSegregated` physics.

## Boundary conditions

[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/energy_inlet_types)
needs one entry for each inlet boundary of the coupled flow, and
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/energy_wall_types)
one entry for each wall boundary, in the same order. Each needs a matching entry in
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/energy_inlet_functors)
or
[!param](/Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics/energy_wall_functors):
a temperature for `fixed-temperature`, or a heat flux into the domain for `heatflux`. Outlets of
the flow need no input: fluid leaving or re-entering through them carries the temperature of the
adjacent cell, and no heat is conducted through them.

## Solver settings

The energy system is solved by the `energy_system` of the `ReducedPressurePIMPLE` executioner,
after the momentum and pressure correction of each outer iteration. The temperature is recovered
right after that solve. Because $\rho h$ is dominated by the large constant part of the phase
enthalpies, a relative linear tolerance on its residual means little: set `energy_l_tol = 0` and
control the linear solve with `energy_l_abs_tol`.

## Example Input Syntax

In this example, a layer of water under a layer of silicone oil sits at rest in a closed tub
between a heated floor and a cooled lid. The water is the liquid phase and the oil the gas phase;
each takes its properties from its own fluid properties object, and both expand with temperature
about the reference temperature `T0`, which defines the buoyancy density `rho_buoyancy`. The four
walls of the flow are, in order, the fixed-temperature floor and lid and the two adiabatic sides.

!listing modules/navier_stokes/test/tests/finite_volume/two_phase/diffuse_interface/oil_water_tub.i block=Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated

The flow physics passes `rho_buoyancy` to the gravity force:

!listing modules/navier_stokes/test/tests/finite_volume/two_phase/diffuse_interface/oil_water_tub.i block=Physics/NavierStokes/ConservativeDiffuseInterfaceFlowSegregated

!syntax parameters /Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics

!syntax inputs /Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics

!syntax children /Physics/NavierStokes/ConservativeDiffuseInterfaceHeatTransferSegregated/WCNSLinearFVConservativeDiffuseInterfaceHeatTransferPhysics
