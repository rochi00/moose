# Analytical Variable-Density Two-Phase Stefan Solution

## Validation solution used by Thirumalaisamy and Bhalla (2023)

This document describes the one-dimensional analytical Stefan solution used to validate the low-Mach enthalpy method in:

> R. Thirumalaisamy and A. P. S. Bhalla, *A low Mach enthalpy method to model non-isothermal gas–liquid–solid flows with melting and solidification*, International Journal of Multiphase Flow 169 (2023) 104605.

The solution differs from the standard Stefan benchmark in three important ways:

1. the solid and liquid densities may be different;
2. the density change generates a uniform liquid velocity through mass conservation; and
3. the interfacial energy balance may retain the kinetic-energy jump associated with the velocity discontinuity.

The analytical fields are the solid and liquid temperatures, the moving interface position, the phase velocities, and the pressure in each phase.

---

## 1. Physical problem

Consider a planar, one-dimensional domain

\[
\Omega = \{x:0\le x\le l\}.
\]

At time \(t=0\), the domain contains only liquid at a uniform temperature

\[
T^L(x,0)=T_i, \qquad T_i>T_m,
\]

where \(T_m\) is the solidification temperature. At \(t=0^+\), the left boundary is suddenly cooled and held at

\[
T^S(0,t)=T_o, \qquad T_o<T_m.
\]

Solid begins to form at the left boundary. The planar solid–liquid interface is located at

\[
x^*=s(t),
\]

and moves in the positive \(x\)-direction. The instantaneous phase domains are

\[
\Omega^S(t)=\{x:0\le x<s(t)\},
\qquad
\Omega^L(t)=\{x:s(t)<x\le l\}.
\]

The analytical thermal boundary condition at the far liquid boundary is

\[
T^L(l,t)=T_i.
\]

The right boundary is also used to set the pressure reference,

\[
p^L(l,t)=p_i,
\]

with \(p_i=0\) in the validation configuration shown schematically in the paper. The domain must be long enough that the liquid temperature has reached the nearly uniform far-field plateau before \(x=l\).

The interface temperature is continuous and fixed at the phase-change temperature:

\[
T^S(s(t)^-,t)=T^L(s(t)^+,t)=T_m.
\]

The interface is planar, so curvature and surface-tension contributions are absent from this benchmark.

---

## 2. Material properties and notation

The phase properties are constant within each bulk phase:

| Quantity | Solid | Liquid |
|---|---:|---:|
| Density | \(\rho^S\) | \(\rho^L\) |
| Specific heat | \(C^S\) | \(C^L\) |
| Thermal conductivity | \(\kappa^S\) | \(\kappa^L\) |
| Thermal diffusivity | \(\alpha^S=\kappa^S/(\rho^S C^S)\) | \(\alpha^L=\kappa^L/(\rho^L C^L)\) |

Define the solid-to-liquid density ratio

\[
R_\rho=\frac{\rho^S}{\rho^L}.
\]

Its physical interpretation is:

- \(R_\rho=1\): no volume change and no density-induced liquid flow;
- \(R_\rho<1\): the solid is less dense, so the material expands during solidification and pushes liquid in the direction of interface motion;
- \(R_\rho>1\): the solid is denser, so the material shrinks during solidification and pulls liquid toward the interface.

The interface speed is

\[
u^*=\dot{s}=\frac{ds}{dt}.
\]

Superscripts \(S\) and \(L\) denote the solid and liquid phases. The notation \(s^-\) and \(s^+\) denotes the limiting locations immediately inside the solid and liquid, respectively.

---

## 3. Density-induced velocity from the mass jump

The Rankine–Hugoniot mass jump at the moving interface is

\[
(\rho^L-\rho^S)u^*
=
\rho^L u^L(s^+,t)-\rho^S u^S(s^-,t).
\]

The solid is assumed stationary:

\[
u^S(x,t)=0.
\]

The liquid is incompressible and one-dimensional, so

\[
\frac{\partial u^L}{\partial x}=0.
\]

Consequently, the liquid velocity is spatially uniform and is determined entirely by the interface speed:

\[
\boxed{
 u^L(x,t)=u^L(s^+,t)=(1-R_\rho)\dot{s}(t)
}
\qquad x\in\Omega^L(t).
\]

Thus, the sharp-interface velocity field is

\[
u(x,t)=
\begin{cases}
0, & 0\le x<s(t),\\[4pt]
(1-R_\rho)\dot{s}(t), & s(t)<x\le l.
\end{cases}
\]

This relation is one of the main differences from the conventional equal-density Stefan problem.

---

## 4. Governing temperature equations

Because the solid is stationary, its temperature satisfies a pure heat equation:

\[
\rho^S C^S\frac{\partial T^S}{\partial t}
=
\kappa^S\frac{\partial^2T^S}{\partial x^2},
\qquad x\in\Omega^S(t).
\]

The liquid temperature satisfies an advection–diffusion equation:

\[
\rho^L C^L
\left(
\frac{\partial T^L}{\partial t}
+u^L\frac{\partial T^L}{\partial x}
\right)
=
\kappa^L\frac{\partial^2T^L}{\partial x^2},
\qquad x\in\Omega^L(t),
\]

where

\[
u^L=(1-R_\rho)\dot{s}.
\]

Equivalently, using thermal diffusivities,

\[
\frac{\partial T^S}{\partial t}
=
\alpha^S\frac{\partial^2T^S}{\partial x^2},
\]

\[
\frac{\partial T^L}{\partial t}
+(1-R_\rho)\dot{s}\frac{\partial T^L}{\partial x}
=
\alpha^L\frac{\partial^2T^L}{\partial x^2}.
\]

---

## 5. Interfacial Stefan condition

### 5.1 Effective latent heat

The energy jump derived in the paper contains latent heat, a specific-heat jump, and a kinetic-energy jump. Define

\[
\boxed{
L_{\mathrm{eff}}
=
L+(C^L-C^S)(T_m-T_r)
}
\]

where

- \(L\) is the latent heat of melting/solidification;
- \(T_r\) is the bulk phase-change temperature measured at the reference pressure used to define phase enthalpies.

For the validation cases in the paper, \(T_r=T_m\), so

\[
L_{\mathrm{eff}}=L.
\]

### 5.2 Full sharp-interface energy balance

The Stefan condition is

\[
\boxed{
\rho^S
\left[
L_{\mathrm{eff}}
-
\frac{1}{2}(1-R_\rho^2)\dot{s}^{\,2}
\right]
\dot{s}
=
\left.
\left(
\kappa^S\frac{\partial T^S}{\partial x}
-
\kappa^L\frac{\partial T^L}{\partial x}
\right)
\right|_{x=s(t)}
}
\]

The term

\[
-\frac{1}{2}(1-R_\rho^2)\dot{s}^{\,2}
\]

is the kinetic-energy contribution per unit mass. It vanishes when \(R_\rho=1\). The full left side contains a term proportional to \(\dot{s}^3\), which becomes potentially important near the initial-time singularity \(\dot{s}\sim t^{-1/2}\).

A Stefan condition derived only from an enthalpy equation generally omits this kinetic-energy jump. The paper derives the condition from total energy precisely to retain it.

---

## 6. Similarity representation

Introduce the dimensionless interface parameter \(\lambda(t)\) through

\[
\boxed{
s(t)=2\lambda(t)\sqrt{\alpha^L t}
}
\]

or, equivalently,

\[
\lambda(t)=\frac{s(t)}{2\sqrt{\alpha^L t}}.
\]

Unlike the standard Stefan parameter, \(\lambda\) can depend on time when the kinetic-energy jump is retained.

### 6.1 Solid similarity variable

Define

\[
\eta^S=\frac{x}{2\sqrt{\alpha^S t}}.
\]

Substitution into the solid heat equation gives

\[
\frac{d^2T^S}{d(\eta^S)^2}
+2\eta^S\frac{dT^S}{d\eta^S}=0,
\]

whose bounded solution is an affine function of \(\operatorname{erf}(\eta^S)\).

### 6.2 Liquid similarity variable

The liquid equation includes the velocity generated by phase-change volume variation. Define

\[
\eta^L
=
\frac{x}{2\sqrt{\alpha^L t}}
-\lambda(t)(1-R_\rho).
\]

This shift follows from choosing the auxiliary function in Appendix A of the paper as

\[
b(t)=-\lambda(t)(1-R_\rho).
\]

The advective contribution then cancels in the similarity transformation, and the liquid equation reduces to

\[
\frac{d^2T^L}{d(\eta^L)^2}
+2\eta^L\frac{dT^L}{d\eta^L}=0.
\]

At the interface \(x=s(t)\), the two similarity variables are

\[
\eta^S_s
=
\lambda(t)\sqrt{\frac{\alpha^L}{\alpha^S}},
\qquad
\eta^L_s
=
\lambda(t)R_\rho.
\]

The factor \(R_\rho\) in the liquid interfacial coordinate is a direct consequence of density-induced advection.

---

## 7. Closed-form temperature fields

The interfacial temperature condition determines the amplitude functions. They are

\[
A(\lambda(t))
=
\frac{T_m-T_o}
{\operatorname{erf}\!\left(
\lambda(t)\sqrt{\alpha^L/\alpha^S}
\right)},
\]

\[
B(\lambda(t))
=
\frac{T_m-T_i}
{\operatorname{erfc}\!\left(\lambda(t)R_\rho\right)}.
\]

The paper emphasizes that \(A\) and \(B\) need not be assumed constant in advance. They may depend on \(\lambda(t)\), and hence implicitly on time, without depending on the similarity coordinates themselves.

### 7.1 Solid temperature

For \(0\le x\le s(t)\),

\[
\boxed{
T^S(x,t)
=
T_o
+(T_m-T_o)
\frac{
\operatorname{erf}\!\left(
\dfrac{x}{2\sqrt{\alpha^S t}}
\right)
}
{
\operatorname{erf}\!\left(
\lambda(t)\sqrt{\dfrac{\alpha^L}{\alpha^S}}
\right)
}
}
\]

This expression satisfies

\[
T^S(0,t)=T_o,
\qquad
T^S(s(t),t)=T_m.
\]

### 7.2 Liquid temperature

For \(s(t)\le x\le l\),

\[
\boxed{
T^L(x,t)
=
T_i
+(T_m-T_i)
\frac{
\operatorname{erfc}\!\left[
\dfrac{x}{2\sqrt{\alpha^L t}}
-
\lambda(t)(1-R_\rho)
\right]
}
{
\operatorname{erfc}\!\left(\lambda(t)R_\rho\right)
}
}
\]

At the interface,

\[
\frac{s(t)}{2\sqrt{\alpha^L t}}
-
\lambda(t)(1-R_\rho)
=
\lambda(t)R_\rho,
\]

so \(T^L(s(t),t)=T_m\). Far ahead of the interface, the complementary error function tends to zero and \(T^L\to T_i\).

---

## 8. Interfacial heat gradients

The solid-side temperature gradient at the interface is

\[
\left.
\frac{\partial T^S}{\partial x}
\right|_{s^-}
=
\frac{T_m-T_o}
{\sqrt{\pi\alpha^S t}\,
\operatorname{erf}\!\left(
\lambda\sqrt{\alpha^L/\alpha^S}
\right)}
\exp\!\left(-\lambda^2\frac{\alpha^L}{\alpha^S}\right).
\]

The liquid-side temperature gradient is

\[
\left.
\frac{\partial T^L}{\partial x}
\right|_{s^+}
=
-
\frac{T_m-T_i}
{\sqrt{\pi\alpha^L t}\,
\operatorname{erfc}(\lambda R_\rho)}
\exp(-\lambda^2R_\rho^2).
\]

These expressions are substituted into the Stefan condition to determine \(\lambda(t)\).

---

## 9. Equation determining \(\lambda(t)\)

Differentiating

\[
s(t)=2\lambda(t)\sqrt{\alpha^L t}
\]

gives the exact interface speed

\[
\boxed{
\dot{s}
=
\lambda\sqrt{\frac{\alpha^L}{t}}
+
2\sqrt{\alpha^L t}\,\dot{\lambda}
}
\]

The full substitution into the Stefan condition produces a differential relation involving time derivatives of \(\lambda\). To retain a convenient analytical evaluation procedure, the paper neglects the \(\dot{\lambda}\)-dependent terms and uses the leading contribution

\[
\dot{s}
\approx
\lambda\sqrt{\frac{\alpha^L}{t}}.
\]

The resulting time-dependent transcendental equation is

\[
\boxed{
\begin{aligned}
&\rho^S
\left[
L_{\mathrm{eff}}
-
\frac{1-R_\rho^2}{2}
\left(
\frac{\lambda^2\alpha^L}{t}
\right)
\right]
\lambda\sqrt{\alpha^L}
\\[4pt]
&\quad=
\kappa^S
\frac{T_m-T_o}
{\operatorname{erf}\!\left(
\lambda\sqrt{\alpha^L/\alpha^S}
\right)}
\frac{
\exp\!\left(-\lambda^2\alpha^L/\alpha^S\right)
}
{\sqrt{\pi\alpha^S}}
\\[4pt]
&\qquad\quad+
\kappa^L
\frac{T_m-T_i}
{\operatorname{erfc}(\lambda R_\rho)}
\frac{\exp(-\lambda^2R_\rho^2)}
{\sqrt{\pi\alpha^L}}.
\end{aligned}
}
\]

Define \(F(\lambda;t)\) as the left side minus the right side. At every required time \(t>0\), the analytical solution is obtained from the physically relevant positive root

\[
F(\lambda;t)=0.
\]

### 9.1 Why \(\lambda\) is time dependent

The kinetic-energy term contains

\[
\frac{\lambda^2\alpha^L}{t},
\]

so retaining it makes the root depend explicitly on time. If the kinetic-energy jump is dropped, the equation becomes

\[
\begin{aligned}
\rho^S L_{\mathrm{eff}}\lambda\sqrt{\alpha^L}
&=
\kappa^S
\frac{T_m-T_o}
{\operatorname{erf}\!\left(
\lambda\sqrt{\alpha^L/\alpha^S}
\right)}
\frac{
\exp\!\left(-\lambda^2\alpha^L/\alpha^S\right)
}
{\sqrt{\pi\alpha^S}}
\\
&\quad+
\kappa^L
\frac{T_m-T_i}
{\operatorname{erfc}(\lambda R_\rho)}
\frac{\exp(-\lambda^2R_\rho^2)}
{\sqrt{\pi\alpha^L}},
\end{aligned}
\]

which is time independent and gives a constant \(\lambda\).

Dropping the kinetic-energy term does **not** by itself impose equal densities: \(R_\rho\) still appears in the liquid temperature and heat-flux terms. The conventional equal-density Stefan solution is recovered only when

\[
R_\rho=1,
\qquad
L_{\mathrm{eff}}=L.
\]

### 9.2 Accuracy of omitting \(\dot{\lambda}\)

For the material properties and density ratios used in the paper, the authors solve the transcendental equation at successive times and then evaluate \(\dot{\lambda}\) in post-processing. They report that

\[
2\sqrt{\alpha^L t}\,\dot{\lambda}
\]

is at least six orders of magnitude smaller than

\[
\lambda\sqrt{\frac{\alpha^L}{t}}.
\]

This supports use of the simpler transcendental equation for the reported validation cases. The approximation is parameter-dependent and should be rechecked for substantially different materials, time scales, or geometries.

---

## 10. Interface position, interface speed, and phase velocities

Once \(\lambda(t)\) is known,

\[
\boxed{
s(t)=2\lambda(t)\sqrt{\alpha^L t}
}
\]

and the exact interface speed is

\[
\dot{s}
=
\lambda\sqrt{\frac{\alpha^L}{t}}
+2\sqrt{\alpha^L t}\,\dot{\lambda}.
\]

The approximation used with the paper's transcendental equation is

\[
\boxed{
\dot{s}\approx\lambda\sqrt{\frac{\alpha^L}{t}}
}
\]

and therefore

\[
\boxed{
u^S=0}
\]

and

\[
\boxed{
u^L=(1-R_\rho)\dot{s}}
\]

throughout the liquid phase.

The sign of \(u^L\) is an immediate diagnostic:

\[
\begin{cases}
u^L>0, & R_\rho<1 \quad\text{(expansion)},\\
u^L=0, & R_\rho=1,\\
u^L<0, & R_\rho>1 \quad\text{(shrinkage)}.
\end{cases}
\]

Both \(\dot{s}\) and \(u^L\) behave as \(t^{-1/2}\) near \(t=0^+\). The analytical solution is therefore evaluated only for \(t>0\); in a discrete calculation, the first evaluation is normally made at \(t=\Delta t\).

---

## 11. Analytical pressure field

### 11.1 Liquid pressure

The one-dimensional liquid momentum equation is

\[
\rho^L
\left(
\frac{\partial u^L}{\partial t}
+u^L\frac{\partial u^L}{\partial x}
\right)
=
-\frac{\partial p^L}{\partial x}
+\mu^L\frac{\partial^2u^L}{\partial x^2}.
\]

Because \(u^L\) is uniform in space,

\[
\frac{\partial u^L}{\partial x}=0,
\qquad
\frac{\partial^2u^L}{\partial x^2}=0,
\]

and therefore

\[
\rho^L(1-R_\rho)\ddot{s}
=
-\frac{\partial p^L}{\partial x}.
\]

The general pressure expression obtained by integrating from \(x\) to the pressure-reference boundary \(l\) is

\[
\boxed{
p^L(x,t)
=
p_i+ho^L(1-R_\rho)\ddot{s}(t)(l-x)
}
\]

which is linear in \(x\).

Using the same leading-order approximation that neglects derivatives of \(\lambda\),

\[
\ddot{s}
\approx
-\frac{\lambda\sqrt{\alpha^L}}{2t^{3/2}},
\]

and the paper writes

\[
\boxed{
p^L(x,t)
=
p_i
-
\frac{\lambda(t)}{2t^{3/2}}
(l-x)(\rho^L-\rho^S)\sqrt{\alpha^L}
}
\]

for \(s(t)<x\le l\).

### 11.2 Pressure jump and solid pressure

Using the jump convention in the paper,

\[
\boxed{
[[p]]
=
p^L(s^+,t)-p^S(s^-,t)
=
\rho^S(1-R_\rho)\dot{s}^{\,2}
}
\]

and the solid pressure is uniform:

\[
\boxed{
p^S(x,t)
=
p^L(s^+,t)-[[p]]
}
\]

or, after substituting the approximate liquid pressure,

\[
\boxed{
p^S(x,t)
=
p_i
-
\frac{\lambda(t)}{2t^{3/2}}
[l-s(t)](\rho^L-\rho^S)\sqrt{\alpha^L}
-[[p]]
}
\]

for \(0\le x<s(t)\).

### 11.3 Important pressure-validation caveat

The analytical pressure solution is a sharp-interface result. The numerical enthalpy method in the paper uses a diffuse mushy region and a Carman–Kozeny drag force to suppress velocity in the solid. That drag produces a much larger pressure variation across the mushy zone than the sharp jump above. Consequently:

- the analytical and numerical solutions should both show a linear liquid pressure and an approximately uniform solid pressure;
- their pressure magnitudes and interfacial jumps are not expected to agree quantitatively for the diffuse-interface formulation used in the paper.

The interface position, temperature, and liquid velocity are the primary quantitative validation fields.

---

## 12. Evaluation algorithm

The following procedure evaluates the analytical solution at a requested time \(t>0\).

### Inputs

Supply

\[
\rho^S,\rho^L,C^S,C^L,\kappa^S,\kappa^L,
L,T_o,T_i,T_m,T_r,l,p_i.
\]

### Derived constants

Compute

\[
R_\rho=\frac{\rho^S}{\rho^L},
\qquad
\alpha^S=\frac{\kappa^S}{\rho^S C^S},
\qquad
\alpha^L=\frac{\kappa^L}{\rho^L C^L},
\]

\[
L_{\mathrm{eff}}
=L+(C^L-C^S)(T_m-T_r).
\]

### Time-by-time procedure

1. **Choose a positive time.** Start at \(t=\Delta t\), not at \(t=0\).
2. **Solve the transcendental equation.** Find the physically continuous positive root \(\lambda(t)\) of \(F(\lambda;t)=0\).
3. **Compute the interface position.**
   \[
   s(t)=2\lambda(t)\sqrt{\alpha^L t}.
   \]
4. **Compute the interface speed.** Use either
   \[
   \dot{s}\approx\lambda\sqrt{\alpha^L/t}
   \]
   or, after obtaining \(\lambda\) at several times, include a finite-difference approximation to \(\dot{\lambda}\).
5. **Compute velocities.**
   \[
   u^S=0,
   \qquad
   u^L=(1-R_\rho)\dot{s}.
   \]
6. **Evaluate temperature piecewise.** Use \(T^S\) for \(x<s(t)\), \(T_m\) at the interface, and \(T^L\) for \(x>s(t)\).
7. **Evaluate pressure if required.** Use the linear liquid pressure, the sharp pressure jump, and the uniform solid pressure.

### Pseudocode

```text
compute R_rho, alpha_S, alpha_L, L_eff

for each requested time t > 0:
    define F(lambda; t) from the transcendental equation
    solve F(lambda; t) = 0 for the positive physical root

    s = 2 * lambda * sqrt(alpha_L * t)

    if lambda_dot is neglected:
        s_dot = lambda * sqrt(alpha_L / t)
    else:
        estimate lambda_dot from neighboring times
        s_dot = lambda * sqrt(alpha_L / t) \
                + 2 * sqrt(alpha_L * t) * lambda_dot

    u_S = 0
    u_L = (1 - R_rho) * s_dot

    for each spatial point x:
        if x < s:
            evaluate T_S(x,t)
        else if x > s:
            evaluate T_L(x,t)
        else:
            T = T_m

    optionally evaluate p_L(x,t), pressure jump, and p_S(x,t)
```

### Root-finding notes

The paper uses MATLAB's `fzero` and solves the equation beginning at \(t=\Delta t\). For a time sequence, the previous root is an effective initial guess for the next time. A bracketed method is preferable when possible.

For numerical stability, the liquid heat-flux ratio may be evaluated as

\[
\frac{e^{-z^2}}{\operatorname{erfc}(z)}
=
\frac{1}{\operatorname{erfcx}(z)},
\qquad z=\lambda R_\rho,
\]

where \(\operatorname{erfcx}\) is the scaled complementary error function. This avoids underflow or division by a very small \(\operatorname{erfc}\) for large positive arguments.

---

## 13. Material values used in the paper

The thermal properties used for the Stefan validations are:

| Property | Symbol | Value |
|---|---:|---:|
| Solid thermal conductivity | \(\kappa^S\) | \(211\ \mathrm{W/(m\,K)}\) |
| Liquid thermal conductivity | \(\kappa^L\) | \(91\ \mathrm{W/(m\,K)}\) |
| Solid specific heat | \(C^S\) | \(910\ \mathrm{J/(kg\,K)}\) |
| Liquid specific heat | \(C^L\) | \(1042.4\ \mathrm{J/(kg\,K)}\) |
| Solidification temperature | \(T_m\) | \(933.6\ \mathrm{K}\) |
| Reference phase-change temperature | \(T_r\) | \(933.6\ \mathrm{K}\) |
| Liquidus temperature used by the diffuse numerical model | \(T_{\mathrm{liq}}\) | \(938.6\ \mathrm{K}\) |
| Solidus temperature used by the diffuse numerical model | \(T_{\mathrm{sol}}\) | \(928.6\ \mathrm{K}\) |
| Latent heat | \(L\) | \(383840\ \mathrm{J/kg}\) |

The sharp analytical solution uses the single interface temperature \(T_m\). The numerical enthalpy method distributes phase change over

\[
\Delta T=T_{\mathrm{liq}}-T_{\mathrm{sol}}=10\ \mathrm{K}.
\]

The common initial and boundary temperatures are

\[
T_i=973.6\ \mathrm{K},
\qquad
T_o=298.6\ \mathrm{K}.
\]

Because \(T_r=T_m\), the validation uses \(L_{\mathrm{eff}}=L\).

---

## 14. Three density cases used for validation

### 14.1 Matched-density case

\[
\rho^S=\rho^L=2475\ \mathrm{kg/m^3},
\qquad
R_\rho=1.
\]

Consequences:

\[
u^L=0,
\]

and the solution reduces to the standard two-phase Stefan solution when \(L_{\mathrm{eff}}=L\).

The paper uses a uniform time step

\[
\Delta t=10^{-3}\ \mathrm{s}.
\]

### 14.2 Volume-expansion case

\[
\rho^S=500\ \mathrm{kg/m^3},
\qquad
\rho^L=2700\ \mathrm{kg/m^3},
\]

\[
R_\rho=\frac{500}{2700}\approx0.185.
\]

The liquid velocity is positive:

\[
u^L=(1-R_\rho)\dot{s}>0.
\]

The density-induced flow carries heat in the direction of front propagation, and the standard equal-density Stefan solution underpredicts the interface position.

The paper uses

\[
\Delta t=10^{-4}\ \mathrm{s}.
\]

### 14.3 Volume-shrinkage case

\[
\rho^S=2700\ \mathrm{kg/m^3},
\qquad
\rho^L=500\ \mathrm{kg/m^3},
\]

\[
R_\rho=5.4.
\]

The liquid velocity is negative:

\[
u^L=(1-R_\rho)\dot{s}<0.
\]

Hot liquid is drawn toward the interface, and the solidification rate is slightly reduced relative to the matched-density reference.

The paper uses the same time step as the expansion case:

\[
\Delta t=10^{-4}\ \mathrm{s}.
\]

Both variable-density cases are simulated to \(t=10\ \mathrm{s}\), as is the matched-density case.

---

## 15. Numerical validation configuration

The low-Mach enthalpy calculation used to validate against the analytical solution is quasi-one-dimensional:

\[
\Omega=[0,1]\times[0,0.05]\ \mathrm{m^2},
\]

with

\[
N_x\times N_y=1280\times64
\]

for the principal comparisons. The domain is periodic in \(y\), and the expected solution is invariant in that direction.

The thermal conditions are:

- initial liquid temperature: \(T_i=973.6\ \mathrm{K}\);
- fixed left temperature: \(T_o=298.6\ \mathrm{K}\);
- adiabatic right boundary in the CFD calculation.

The analytical solution uses \(T^L(l,t)=T_i\). In the long numerical domain, the right boundary lies in the nearly uniform liquid-temperature plateau, so the adiabatic numerical condition approximates the analytical far-field condition over the validation interval.

The flow boundary conditions are:

- zero velocity at the left boundary;
- zero pressure/outflow at the right boundary.

The solid and liquid viscosities are set to zero for these cases. The numerical solid–liquid interface is extracted from the liquid-fraction contour

\[
\varphi=0.5.
\]

---

## 16. Quantities compared in the paper

The principal analytical-to-CFD comparisons are:

### 16.1 Interface trajectory

\[
x^*(t)=s(t)=2\lambda(t)\sqrt{\alpha^L t}.
\]

This is compared over the full simulation interval \(0<t\le10\ \mathrm{s}\).

### 16.2 Temperature profiles

The piecewise analytical temperature

\[
T(x,t)=
\begin{cases}
T^S(x,t), & x<s(t),\\[3pt]
T_m, & x=s(t),\\[3pt]
T^L(x,t), & x>s(t)
\end{cases}
\]

is compared with numerical profiles at representative times, including \(t=1\), \(5\), and \(10\ \mathrm{s}\).

### 16.3 Liquid velocity

The analytical velocity is uniform in the liquid:

\[
u^L(t)=(1-R_\rho)\dot{s}(t).
\]

The matched-density case has no flow. Expansion produces positive velocity, while shrinkage produces negative velocity.

### 16.4 Pressure structure

The sharp analytical result predicts

- linear pressure in the liquid;
- uniform pressure in the solid;
- a discontinuous interfacial pressure jump.

Only the spatial structure is directly comparable to the diffuse enthalpy calculation because the Carman–Kozeny mushy-zone pressure drop changes the numerical magnitude.

### 16.5 Error measure used in the convergence study

For a sampled quantity \(\psi\), the paper defines a root-mean-square error

\[
\|\mathcal{E}_\psi\|_{\mathrm{RMSE}}
=
\frac{
\|\psi_{\mathrm{reference}}-\psi_{\mathrm{numerical}}\|_2
}
{\sqrt{\mathcal{N}}},
\]

where \(\mathcal{N}\) is the number of samples. The analytical solution is used as one reference; the finest-grid diffuse solution is used as a second reference to distinguish sharp-interface error from the discretization error of the diffuse model itself.

---

## 17. Reference values computed from the reported parameters

The following values are obtained by directly evaluating the equations above with the material data reported in the paper. They are provided as implementation sanity checks; they are computed values rather than a table reproduced from the article.

### 17.1 Thermal diffusivities and \(\lambda\)

| Case | \(R_\rho\) | \(\alpha^S\;[\mathrm{m^2/s}]\) | \(\alpha^L\;[\mathrm{m^2/s}]\) | \(\lambda\) at \(t=10\,\mathrm{s}\) |
|---|---:|---:|---:|---:|
| Matched density | \(1\) | \(9.368409\times10^{-5}\) | \(3.527214\times10^{-5}\) | \(1.12533235\) |
| Expansion | \(0.185185\) | \(4.637363\times10^{-4}\) | \(3.233279\times10^{-5}\) | \(2.55867445\) |
| Shrinkage | \(5.4\) | \(8.587709\times10^{-5}\) | \(1.745971\times10^{-4}\) | \(0.48781783\) |

For these parameters, the kinetic-energy correction changes \(\lambda\) only slightly over the reported time interval.

### 17.2 Interface and velocity at \(t=10\,\mathrm{s}\)

Using \(\dot{s}\approx\lambda\sqrt{\alpha^L/t}\):

| Case | \(s(10\,\mathrm{s})\;[\mathrm{m}]\) | \(\dot{s}(10\,\mathrm{s})\;[\mathrm{m/s}]\) | \(u^L(10\,\mathrm{s})\;[\mathrm{m/s}]\) |
|---|---:|---:|---:|
| Matched density | \(0.0422695\) | \(0.00211347\) | \(0\) |
| Expansion | \(0.0920167\) | \(0.00460083\) | \(0.00374883\) |
| Shrinkage | \(0.0407668\) | \(0.00203834\) | \(-0.00896869\) |

These values reproduce the qualitative ordering shown in the paper: expansion advances the front substantially, while shrinkage gives a front position close to but slightly below the matched-density solution.

---

## 18. Important modeling and implementation cautions

### 18.1 The solution is defined for \(t>0\)

At \(t=0\), no solid exists, \(T^S\) is not defined, and

\[
s(t)=2\lambda(t)\sqrt{\alpha^L t}
\]

has a singular derivative. Begin analytical and numerical comparisons at a positive time.

### 18.2 The analytical interface is sharp

The analytical problem has a single temperature \(T_m\) at an infinitesimally thin interface. The enthalpy calculation uses a finite mushy interval \(T_{\mathrm{sol}}\le T\le T_{\mathrm{liq}}\). Agreement should improve as the numerical mushy interval and grid spacing are reduced together, while still resolving the mushy zone.

### 18.3 Do not omit liquid advection for unequal densities

For \(R_\rho\ne1\), using a diffusion-only liquid heat equation is inconsistent with the interfacial mass jump. The liquid advection term is essential and is responsible for the shifted liquid similarity coordinate.

### 18.4 Distinguish two simplifications

Two different reductions are often conflated:

1. **Omit the kinetic-energy jump:** \(\lambda\) becomes constant, but density-induced liquid advection remains if \(R_\rho\ne1\).
2. **Assume equal densities:** \(R_\rho=1\), so liquid velocity and all density-change effects vanish.

### 18.5 Check the \(\dot{\lambda}\) approximation for new parameters

The paper verifies that \(\dot{\lambda}\) is negligible for its validation cases. That observation should not automatically be transferred to dramatically different latent heats, density ratios, length scales, or very early times.

### 18.6 Pressure is formulation-sensitive

The sharp pressure jump is not expected to match the pressure drop across a diffuse, drag-regularized mushy layer. Use pressure primarily to verify the expected phasewise spatial form unless the numerical method imposes the sharp jump directly.

---

## 19. Compact formula sheet

### Definitions

\[
R_\rho=\frac{\rho^S}{\rho^L},
\qquad
\alpha^S=\frac{\kappa^S}{\rho^S C^S},
\qquad
\alpha^L=\frac{\kappa^L}{\rho^L C^L},
\]

\[
L_{\mathrm{eff}}=L+(C^L-C^S)(T_m-T_r).
\]

### Interface

\[
s=2\lambda\sqrt{\alpha^L t},
\]

\[
\dot{s}
=
\lambda\sqrt{\frac{\alpha^L}{t}}
+2\sqrt{\alpha^L t}\,\dot{\lambda}
\approx
\lambda\sqrt{\frac{\alpha^L}{t}}.
\]

### Velocities

\[
u^S=0,
\qquad
u^L=(1-R_\rho)\dot{s}.
\]

### Temperatures

\[
T^S
=
T_o+(T_m-T_o)
\frac{
\operatorname{erf}\left(x/(2\sqrt{\alpha^S t})\right)
}
{
\operatorname{erf}\left(\lambda\sqrt{\alpha^L/\alpha^S}\right)
},
\]

\[
T^L
=
T_i+(T_m-T_i)
\frac{
\operatorname{erfc}\left(x/(2\sqrt{\alpha^L t})-\lambda(1-R_\rho)\right)
}
{
\operatorname{erfc}(\lambda R_\rho)
}.
\]

### Stefan condition

\[
\rho^S
\left[
L_{\mathrm{eff}}
-
\frac{1}{2}(1-R_\rho^2)\dot{s}^{\,2}
\right]
\dot{s}
=
\left.
\left(
\kappa^S T_x^S-\kappa^L T_x^L
\right)
\right|_{x=s}.
\]

### Pressure

\[
p^L=p_i+\rho^L(1-R_\rho)\ddot{s}(l-x),
\]

\[
[[p]]=p^L-p^S=\rho^S(1-R_\rho)\dot{s}^{\,2},
\]

\[
p^S=p^L(s^+,t)-[[p]].
\]

---

## 20. Source locations in the paper

- Jump conditions and the kinetic-energy-corrected Stefan condition: Section 2, especially Eqs. (4)–(12).
- Variable-density Stefan problem and closed-form solution: Section 3, Eqs. (13)–(29).
- Similarity-variable derivation: Appendix A, Eqs. (A.1)–(A.11).
- Material data and validation configuration: Section 6.1 and Table 1.
- Comparisons for matched density, expansion, and shrinkage: Sections 6.1.1–6.1.3 and Fig. 2.
- Assessment of the neglected \(\dot{\lambda}\) term: Section 6.1.6 and Fig. 8.
- Pressure comparison and diffuse-interface caveat: Section 6.1.7 and Fig. 9.

