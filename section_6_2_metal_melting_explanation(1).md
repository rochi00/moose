# Explanation of Sections 6.2 and 6.3: Metal Melting and Solidification with a Free Surface

This document explains the two three-phase validation and sanity-check cases presented in Sections 6.2 and 6.3 of the paper. Section 6.2 treats melting-induced contraction, while Section 6.3 treats solidification-induced shrinkage and expansion defects.

---

## Section 6.2: Metal Melting with a Free Surface

### 1. Purpose of the test

Section 6.2 presents a **three-phase metal-melting benchmark** for the low-Mach enthalpy method developed in the paper. Unlike the one-dimensional Stefan problems in Section 6.1, this test is not compared with a full analytical transient solution. Instead, it is designed as a controlled numerical sanity check for two specific capabilities:

1. capturing the volume change of a phase-change material in the presence of a gas phase; and
2. handling the complete disappearance of one material phase from the computational domain.

The calculation tests whether the method can correctly represent the following sequence:

\[
\text{solid metal melts}
\;\longrightarrow\;
\text{density changes}
\;\longrightarrow\;
\text{metal volume changes}
\;\longrightarrow\;
\text{the gas--metal free surface moves}.
\]

The example is inspired by a related phase-field calculation by Huang et al., but the authors use their low-Mach enthalpy formulation together with aluminum-based thermophysical properties.

---

### 2. Computational domain and grid

The computational domain is the unit square

\[
\Omega=[0,1]\times[0,1],
\]

which is discretized using

\[
N_x\times N_y=256\times256
\]

uniform cells.

The configuration is effectively one-dimensional in its initial phase arrangement because the three phases occupy horizontal layers, although the governing equations are solved in two dimensions.

---

### 3. Initial phase configuration

At \(t=0\), the domain contains three horizontal regions:

\[
\begin{cases}
0\le y\le0.30\ \mathrm{m}: & \text{liquid metal},\\[2mm]
0.30<y\le0.45\ \mathrm{m}: & \text{solid metal},\\[2mm]
0.45<y\le1.00\ \mathrm{m}: & \text{gas}.
\end{cases}
\]

Therefore,

- the initial liquid depth is \(0.30\ \mathrm{m}\);
- the initial solid-layer thickness is \(0.15\ \mathrm{m}\); and
- the initial gas--metal interface is located at

\[
y_\Gamma(0)=0.45\ \mathrm{m}.
\]

Figure 10 of the paper shows liquid in blue, solid in yellow, and gas in white.

---

### 4. Material properties

#### 4.1 Liquid and solid densities

The liquid is denser than the solid:

\[
\rho^L=2700\ \mathrm{kg/m^3},
\qquad
\rho^S=2475\ \mathrm{kg/m^3}.
\]

Thus,

\[
\rho^L>\rho^S.
\]

This density ordering is the central physical feature of the test. When a fixed mass of the lighter solid melts into the denser liquid, it occupies a smaller volume. The total metal volume must therefore decrease, causing the gas--metal free surface to move downward.

#### 4.2 Viscosity

The liquid and solid are assigned the same viscosity:

\[
\mu^L=\mu^S=1.4\times10^{-3}\ \mathrm{kg/(m\,s)}.
\]

The paper notes that the solid viscosity is fictitious and does not determine the result. Motion in the solid region is instead suppressed by the Carman--Kozeny drag term in the momentum equation.

#### 4.3 Remaining metal properties

The remaining solid and liquid properties are taken from Table 1 of the paper:

\[
\begin{aligned}
\kappa^S &=211\ \mathrm{W/(m\,K)},\\
\kappa^L &=91\ \mathrm{W/(m\,K)},\\
C^S &=910\ \mathrm{J/(kg\,K)},\\
C^L &=1042.4\ \mathrm{J/(kg\,K)},\\
T_m &=933.6\ \mathrm{K},\\
T_{\mathrm{sol}}&=928.6\ \mathrm{K},\\
T_{\mathrm{liq}}&=938.6\ \mathrm{K},\\
L&=383840\ \mathrm{J/kg}.
\end{aligned}
\]

#### 4.4 Gas properties

The gas is assigned air-based properties:

\[
\begin{aligned}
\rho^G &=0.4\ \mathrm{kg/m^3},\\
\kappa^G &=6.1\times10^{-2}\ \mathrm{W/(m\,K)},\\
C^G &=1100\ \mathrm{J/(kg\,K)},\\
\mu^G &=4\times10^{-5}\ \mathrm{kg/(m\,s)}.
\end{aligned}
\]

The liquid-to-gas density ratio is approximately

\[
\frac{\rho^L}{\rho^G}
=
\frac{2700}{0.4}
\approx6750.
\]

Consequently, this example also tests the robustness of the algorithm for a very large gas--metal density ratio.

---

### 5. Initial and boundary temperatures

The initial temperature is specified separately in each phase:

\[
T(\mathbf{x},0)=
\begin{cases}
5T_m, & \text{liquid},\\
0.9T_m, & \text{solid},\\
0.9T_m, & \text{gas}.
\end{cases}
\]

Using \(T_m=933.6\ \mathrm{K}\), these values are

\[
5T_m=4668\ \mathrm{K},
\qquad
0.9T_m=840.24\ \mathrm{K}.
\]

The bottom wall is maintained at

\[
T=6T_m=5601.6\ \mathrm{K}.
\]

Thus, the initially hot liquid and the heated bottom wall provide the thermal energy needed to melt the solid layer.

The boundary conditions are:

- periodic boundary conditions in the \(x\)-direction;
- zero velocity at the bottom wall;
- fixed temperature \(T=6T_m\) at the bottom wall;
- zero pressure/outflow at the top boundary; and
- zero heat flux at the top boundary.

The paper labels both the top and bottom boundaries as \(y=0\). Since the domain is \([0,1]^2\) and the bottom wall is explicitly located at \(y=0\), the top-boundary label appears to be a typographical inconsistency; geometrically, the top boundary is at \(y=1\).

The very high prescribed temperatures should be understood as part of a numerical melting benchmark. Evaporation and condensation are not modeled in this example.

---

### 6. Representation of the two interfaces

The calculation contains two physically different interfaces, and the method represents them with two different variables.

#### 6.1 Gas--metal interface

The outer gas--PCM interface is represented by the zero contour of a level-set function:

\[
\phi(\mathbf{x},t)=0.
\]

The corresponding smoothed Heaviside function \(H\) distinguishes the gas and PCM regions:

\[
H\approx
\begin{cases}
0, & \text{gas},\\
1, & \text{solid or liquid PCM}.
\end{cases}
\]

The level-set function is advected by the computed velocity:

\[
\frac{\partial\phi}{\partial t}
+
\mathbf{u}\cdot\nabla\phi
=0.
\]

The gas--metal interface therefore moves in response to the velocity generated by phase-change-induced volume contraction.

#### 6.2 Solid--liquid interface

The internal solid--liquid transition is represented by the liquid fraction \(\varphi\):

\[
\varphi=
\begin{cases}
0, & \text{solid},\\
1, & \text{liquid},\\
0<\varphi<1, & \text{mushy region}.
\end{cases}
\]

The nonlinear enthalpy equation determines the temperature and enthalpy fields, from which \(\varphi\) is evaluated. As heat enters the solid layer, \(\varphi\) increases from zero to one until the solid phase disappears.

The roles of the two variables are therefore distinct:

- \(\phi\) tracks the **outer gas--metal free surface**;
- \(\varphi\) tracks the **internal solid--liquid phase change**.

---

### 7. Why melting causes volume contraction

The low-Mach formulation imposes the velocity-divergence constraint

\[
\nabla\cdot\mathbf{u}
=
\frac{\rho^S-\rho^L}{\rho}
H\frac{D\varphi}{Dt}.
\]

For this case,

\[
\rho^S-\rho^L
=
2475-2700
=
-225\ \mathrm{kg/m^3}.
\]

During melting,

\[
\frac{D\varphi}{Dt}>0.
\]

It follows that

\[
\nabla\cdot\mathbf{u}<0
\]

inside the mushy region.

A negative velocity divergence represents local volume contraction. The contraction generates a velocity field that pulls the gas--metal free surface downward.

Once all of the solid has melted,

\[
\varphi=1
\]

throughout the PCM and

\[
\frac{D\varphi}{Dt}=0.
\]

The phase-change volume source then vanishes:

\[
\nabla\cdot\mathbf{u}=0.
\]

The metal volume should therefore stop changing after complete melting. The simulation is deliberately continued beyond the time of phase disappearance to verify that the numerical method satisfies this requirement.

---

### 8. Expected final free-surface position

The most important quantitative prediction follows directly from conservation of metal mass.

Because the domain has unit width, the initial metal mass per unit out-of-plane depth is

\[
m_0
=
\rho^L(0.30)
+
\rho^S(0.15).
\]

Substituting the densities gives

\[
\begin{aligned}
m_0
&=2700(0.30)+2475(0.15)\\
&=810+371.25\\
&=1181.25\ \mathrm{kg/m}.
\end{aligned}
\]

After complete melting, all metal is liquid. Let the final liquid depth be \(h_f\). The final metal mass is

\[
m_f=\rho^L h_f.
\]

Conservation of mass requires

\[
\rho^L h_f=m_0,
\]

so

\[
h_f
=
\frac{1181.25}{2700}
=
0.4375\ \mathrm{m}.
\]

The same result can be obtained by converting the original solid thickness into an equivalent liquid thickness:

\[
0.15\frac{\rho^S}{\rho^L}
=
0.15\frac{2475}{2700}
=
0.1375\ \mathrm{m}.
\]

Adding this to the original liquid depth gives

\[
0.30+0.1375=0.4375\ \mathrm{m}.
\]

Therefore, the free surface should fall from

\[
0.45\ \mathrm{m}
\]

to

\[
\boxed{0.4375\ \mathrm{m}}.
\]

The expected downward displacement is

\[
\Delta h
=
0.45-0.4375
=
0.0125\ \mathrm{m}.
\]

Relative to the initial metal depth, this is

\[
\frac{0.0125}{0.45}\times100
\approx2.78\%.
\]

This final interface position is the principal sanity check in Section 6.2.

---

### 9. Time integration

The simulation uses a uniform time step

\[
\Delta t=10^{-3}\ \mathrm{s}
\]

and is continued until

\[
t=250\ \mathrm{s}.
\]

The four snapshots in Figure 10 correspond to \(t=0\), \(100\), \(150\), and \(250\ \mathrm{s}\).

---

### 10. Evolution shown in Figure 10

#### 10.1 At \(t=0\)

The initial three-layer structure is present:

- liquid below \(y=0.30\);
- solid between \(y=0.30\) and \(y=0.45\); and
- gas above \(y=0.45\).

#### 10.2 At \(t=100\ \mathrm{s}\)

Most of the solid has melted. A thin solid layer remains above the liquid. The free surface has already moved downward because the newly formed liquid is denser than the solid from which it formed.

#### 10.3 At \(t=150\ \mathrm{s}\)

The solid phase has disappeared completely. The PCM is entirely liquid, and the gas--liquid interface is located at approximately

\[
y=0.4375\ \mathrm{m},
\]

in agreement with the mass-conservation prediction.

#### 10.4 At \(t=250\ \mathrm{s}\)

The phase configuration is essentially unchanged from \(t=150\ \mathrm{s}\). There is

- no reappearance of solid;
- no continued recession of the free surface; and
- no spurious interfacial dynamics.

Continuing the calculation for another \(100\ \mathrm{s}\) after complete melting demonstrates that the phase-change volume source turns off when the mushy region disappears.

---

### 11. Metal-mass conservation

The total PCM mass is evaluated as

\[
m(t)
=
\int_\Omega
\left[
\rho^L H\varphi
+
\rho^S\left(H-H\varphi\right)
\right]
\,d\Omega.
\]

The first term represents liquid metal, while the second represents solid metal.

The relative percentage mass error is defined as

\[
\mathcal{E}(t)
=
\frac{|m(t)-m_0|}{m_0}\times100.
\]

For the \(256\times256\) calculation in Section 6.2, the reported metal-mass change is approximately

\[
\mathcal{E}\approx0.027\%.
\]

The paper attributes the nonzero error mainly to two sources:

1. the non-conservative character of the level-set method used to track the gas--PCM interface; and
2. numerical error in evaluating the right-hand side of the low-Mach divergence equation.

Appendix C shows that the mass error decreases under grid refinement. On the finer grids considered, it approaches approximately

\[
0.007\%.
\]

The authors also state that reducing the time-step size lowers the error, although those data are not shown.

---

### 12. What Section 6.2 validates

The test verifies the complete physical and numerical chain

\[
\text{heat transfer}
\longrightarrow
\text{solid melting}
\longrightarrow
\varphi\ \text{increases}
\longrightarrow
\nabla\cdot\mathbf{u}<0
\longrightarrow
\text{metal contracts}
\longrightarrow
\text{the free surface moves downward}.
\]

More specifically, Section 6.2 demonstrates:

- coupling of a solid--liquid PCM to a gas phase;
- conversion of solid material into liquid;
- complete disappearance of the solid phase;
- phase-change-induced volume contraction;
- movement of the free surface by a non-divergence-free velocity field;
- recovery of the expected final free-surface position;
- termination of volume change after melting is complete;
- numerical stability for a metal-to-gas density ratio of several thousand; and
- small, grid-convergent PCM-mass error.

---

### 13. What the example does not validate

Section 6.2 does not provide:

- an analytical transient solution for the melting-front trajectory;
- comparison with an aluminum-melting experiment;
- quantitative errors for temperature profiles;
- a formal convergence rate for the melting-front position; or
- a physical evaporation model for the very high prescribed temperatures.

The example is therefore best interpreted as a **three-phase conservation, volume-change, and phase-disappearance benchmark**, rather than as a comprehensive experimental validation of aluminum melting.

---

### 14. Main conclusion

The central result is that conservation of mass predicts the final gas--liquid interface position to be

\[
\boxed{y_\Gamma^{\mathrm{final}}=0.4375\ \mathrm{m}}.
\]

The numerical calculation reproduces this contraction, completely removes the solid phase by approximately \(t=150\ \mathrm{s}\), and remains stationary through \(t=250\ \mathrm{s}\). This behavior confirms that the low-Mach enthalpy formulation activates volume change only while phase change is occurring and turns the associated divergence source off once the material is fully liquid.

---

## Section 6.3: Metal Solidification and Casting Defects

### 1. Purpose of the test

Section 6.3 considers the reverse process of Section 6.2: initially molten aluminum solidifies while exposed to a gas phase. The example is designed to demonstrate why a phase-change model with different liquid and solid densities must also permit a nonzero velocity divergence in the phase-change region.

The principal question is whether the numerical method can reproduce a common casting defect known as a **shrinkage pipe**. When the solid is denser than the liquid, a fixed mass occupies less volume after solidification. Gas must enter the volume vacated by the contracting metal, causing the free surface to cave into the casting and leaving a pipe-shaped depression.

The authors use this example to compare two formulations:

1. a conventional calculation in which velocity is constrained by

   \[
   \nabla\cdot\mathbf{u}=0,
   \]

   even though the liquid and solid densities differ; and
2. the proposed low-Mach enthalpy calculation, in which the divergence is nonzero only in the mushy phase-change region.

The comparison isolates the role of the low-Mach kinematic constraint. It shows that merely placing variable density in the momentum and energy equations is not sufficient to create the volume change required by mass conservation.

---

### 2. Physical meaning of a shrinkage pipe

For a fixed metal mass,

\[
M=\rho V.
\]

If the density increases during solidification, then the occupied volume must decrease:

\[
\rho^S>\rho^L
\quad\Longrightarrow\quad
V^S=\frac{M}{\rho^S}<\frac{M}{\rho^L}=V^L.
\]

In a casting with an exposed upper surface, the gas phase fills the volume lost by the metal. Because cooling begins from the side walls in this problem, solid develops along both sides while a liquid channel remains near the center. As contraction proceeds, the gas-metal interface is drawn downward over this central region. The final solid therefore contains a narrow, pipe-like surface depression rather than a flat upper boundary.

The pipe is not introduced geometrically or imposed as a boundary condition. It is an emergent consequence of

\[
\text{side-wall cooling}
\longrightarrow
\text{solidification}
\longrightarrow
\text{density increase}
\longrightarrow
\text{negative velocity divergence}
\longrightarrow
\text{free-surface recession}.
\]

---

### 3. Computational domain and grid

The calculation is performed in a square domain

\[
\Omega=[0,8\times10^{-3}]\times[0,8\times10^{-3}]\ \mathrm{m^2}.
\]

Thus, the domain is an \(8\ \mathrm{mm}\) by \(8\ \mathrm{mm}\) square. It is discretized using

\[
N_x\times N_y=256\times256
\]

uniform cells. The corresponding cell width is

\[
\Delta x=\Delta y
=\frac{8\times10^{-3}}{256}
=3.125\times10^{-5}\ \mathrm{m}
=31.25\ \mu\mathrm{m}.
\]

The paper identifies the corresponding IBAMR benchmark as

```text
examples/phase_change/ex4
```

---

### 4. Initial phase arrangement

Initially, there is no solid metal in the domain. The phase arrangement is

\[
\begin{cases}
0\le y\le5\times10^{-3}\ \mathrm{m}: & \text{liquid aluminum},\\[2mm]
5\times10^{-3}<y\le8\times10^{-3}\ \mathrm{m}: & \text{gas}.
\end{cases}
\]

Therefore,

- the initial liquid depth is \(5\ \mathrm{mm}\);
- the initial gas depth is \(3\ \mathrm{mm}\); and
- the initially flat gas-liquid interface is at

\[
y_\Gamma(0)=5\times10^{-3}\ \mathrm{m}.
\]

Figure 11(A), Figure 12(A), and Figure 13(A) all begin from this same phase distribution: blue liquid below and white gas above.

---

### 5. Thermophysical properties

#### 5.1 Physical shrinkage case

For the aluminum solidification case that forms a pipe defect, the densities are

\[
\rho^L=2475\ \mathrm{kg/m^3},
\qquad
\rho^S=2700\ \mathrm{kg/m^3}.
\]

Hence,

\[
\rho^S>\rho^L,
\]

so the metal contracts during solidification.

The other solid, liquid, and gas properties are taken from Section 6.2 and Table 1:

\[
\begin{aligned}
\kappa^S &=211\ \mathrm{W/(m\,K)},
&\kappa^L &=91\ \mathrm{W/(m\,K)},\\
C^S &=910\ \mathrm{J/(kg\,K)},
&C^L &=1042.4\ \mathrm{J/(kg\,K)},\\
T_m &=933.6\ \mathrm{K},
&T_{\mathrm{sol}}&=928.6\ \mathrm{K},\\
T_{\mathrm{liq}}&=938.6\ \mathrm{K},
&L&=383840\ \mathrm{J/kg},\\
\mu^L&=\mu^S=1.4\times10^{-3}\ \mathrm{kg/(m\,s)}.
\end{aligned}
\]

The gas properties are

\[
\begin{aligned}
\rho^G &=0.4\ \mathrm{kg/m^3},\\
\kappa^G &=6.1\times10^{-2}\ \mathrm{W/(m\,K)},\\
C^G &=1100\ \mathrm{J/(kg\,K)},\\
\mu^G &=4\times10^{-5}\ \mathrm{kg/(m\,s)}.
\end{aligned}
\]

The liquid-aluminum/gas surface-tension coefficient is

\[
\sigma=0.87\ \mathrm{N/m}.
\]

The momentum equation also includes variable density and gravity. Section 6.3 does not separately state a numerical value for the gravitational acceleration.

#### 5.2 Hypothetical expansion case

The authors also reverse the liquid and solid densities:

\[
\rho^L=2700\ \mathrm{kg/m^3},
\qquad
\rho^S=2475\ \mathrm{kg/m^3}.
\]

This second density ordering is intentionally hypothetical for the test. It causes the material to expand rather than contract when it solidifies.

---

### 6. Initial temperature and thermal boundary conditions

The initial liquid temperature is

\[
T^L(\mathbf{x},0)=2T_m,
\]

whereas the initial gas temperature is

\[
T^G(\mathbf{x},0)=0.5T_m.
\]

Using \(T_m=933.6\ \mathrm{K}\),

\[
2T_m=1867.2\ \mathrm{K},
\qquad
0.5T_m=466.8\ \mathrm{K}.
\]

The thermal boundary conditions are:

- fixed temperature \(T=0.5T_m\) on the left, right, and top boundaries;
- zero heat flux on the bottom boundary.

Because the side walls are cold, solidification begins at the left and right boundaries and progresses inward. The paper states that the top boundary has little influence on solidification because the gas has low thermal conductivity. The insulated bottom prevents cooling from below, so the phase pattern is dominated by lateral solidification.

---

### 7. Velocity and pressure boundary conditions

The flow boundary conditions are:

- zero pressure/outflow at the top boundary;
- zero velocity on the left, right, and bottom boundaries.

The open top permits the gas region and free surface to respond to metal-volume change, while the remaining walls behave as stationary solid boundaries.

---

### 8. Time integration

The simulation uses a uniform time step

\[
\Delta t=10^{-5}\ \mathrm{s}
\]

and is run until

\[
t=1\ \mathrm{s}.
\]

The figures show phase distributions at

\[
t=0,\ 0.1,\ 0.2,\ \text{and}\ 1\ \mathrm{s}.
\]

The authors report that solidification is essentially complete at approximately

\[
t=0.25\ \mathrm{s}.
\]

The remainder of the simulation tests whether the fully solidified configuration remains stationary without spurious phase change or free-surface motion.

---

### 9. Low-Mach divergence during solidification

The low-Mach enthalpy formulation imposes

\[
\nabla\cdot\mathbf{u}
=
\frac{\rho^S-\rho^L}{\rho}
H\frac{D\varphi}{Dt},
\]

where \(\varphi\) is the liquid fraction. During solidification,

\[
\frac{D\varphi}{Dt}<0.
\]

#### 9.1 Shrinkage case

For the physical casting case,

\[
\rho^S-\rho^L
=2700-2475
=225\ \mathrm{kg/m^3}>0.
\]

Therefore,

\[
\frac{D\varphi}{Dt}<0
\quad\Longrightarrow\quad
\nabla\cdot\mathbf{u}<0.
\]

Negative divergence represents local volume contraction in the mushy region. The resulting velocity field draws the free surface downward and allows gas to occupy the lost metal volume.

#### 9.2 Expansion case

When the densities are reversed,

\[
\rho^S-\rho^L
=2475-2700
=-225\ \mathrm{kg/m^3}<0.
\]

Since \(D\varphi/Dt<0\) during solidification,

\[
\nabla\cdot\mathbf{u}>0.
\]

Positive divergence represents local expansion, so the free surface rises and produces a protrusion.

#### 9.3 After solidification is complete

Once the metal is entirely solid,

\[
\varphi=0,
\qquad
\frac{D\varphi}{Dt}=0,
\]

and therefore

\[
\nabla\cdot\mathbf{u}=0.
\]

The volume-change source automatically switches off. This is why the final pipe or protrusion remains stationary after approximately \(0.25\ \mathrm{s}\).

---

### 10. Divergence-free comparison: Figure 11

The first calculation artificially enforces

\[
\nabla\cdot\mathbf{u}=0
\]

throughout the domain.

The liquid and solid densities are still different, and variable density and gravity remain in the momentum equation. Nevertheless, the metal cannot change its total volume because a globally divergence-free velocity field admits no phase-change volume source.

Figure 11 shows the consequences:

#### At \(t=0\)

The liquid occupies the lower \(5\ \mathrm{mm}\), and the free surface is flat.

#### At \(t=0.1\ \mathrm{s}\)

Solid regions have developed along the cold left and right walls. A broad liquid region remains in the center. The free surface remains flat.

#### At \(t=0.2\ \mathrm{s}\)

The side solidification fronts have moved farther inward, leaving a narrow liquid channel near the center. The upper metal surface still does not deform.

#### At \(t=1\ \mathrm{s}\)

The liquid has solidified completely, but the solid occupies exactly the same gross region as the initial liquid. The gas-solid interface remains at its original height.

This result is unphysical for \(\rho^S>\rho^L\). It demonstrates that including different phase densities in constitutive properties and body forces does not by itself enforce the mass-conserving volume change associated with phase transition.

---

### 11. Non-divergence-free shrinkage calculation: Figure 12

The second calculation uses the low-Mach divergence constraint. In the mushy region, solidification produces \(\nabla\cdot\mathbf{u}<0\), so the metal contracts.

#### At \(t=0\)

The initial gas-liquid interface is flat at \(y=5\ \mathrm{mm}\).

#### At \(t=0.1\ \mathrm{s}\)

Solid has formed along the left and right walls. The free surface has already begun to sag near the middle because the partially solidified metal occupies less volume.

#### At \(t=0.2\ \mathrm{s}\)

Only a narrow liquid channel remains near the center. The gas-metal interface has been pulled sharply downward into this central region. The developing geometry has the characteristic shape of a shrinkage pipe.

#### At \(t=1\ \mathrm{s}\)

The metal is completely solid. A narrow gas-filled depression extends downward into the solidified casting. The pipe remains after phase change ends, and the paper reports no spurious gas-solid interfacial motion beyond approximately \(t=0.25\ \mathrm{s}\).

The contrast between Figures 11 and 12 is the central result of Section 6.3:

\[
\boxed{
\text{different densities} + \nabla\cdot\mathbf{u}=0
\;\not\Rightarrow\;
\text{correct volume change}
}
\]

whereas

\[
\boxed{
\text{different densities} +
\text{low-Mach phase-change divergence}
\;\Rightarrow\;
\text{pipe formation}
}.
\]

---

### 12. Derived global shrinkage estimate

The following calculation is not explicitly written in Section 6.3, but it follows directly from the dimensions and densities reported in the paper and provides a useful global mass-conservation check.

Per unit out-of-plane depth, the initial liquid area is

\[
A_0
=(8\times10^{-3})(5\times10^{-3})
=4.0\times10^{-5}\ \mathrm{m^2}.
\]

The initial metal mass per unit depth is

\[
\begin{aligned}
m_0
&=\rho^L A_0\\
&=2475(4.0\times10^{-5})\\
&=9.90\times10^{-2}\ \mathrm{kg/m}.
\end{aligned}
\]

After complete solidification, conservation of mass requires

\[
A_f=\frac{m_0}{\rho^S}
=A_0\frac{\rho^L}{\rho^S}.
\]

Thus,

\[
A_f
=(4.0\times10^{-5})\frac{2475}{2700}
=3.6667\times10^{-5}\ \mathrm{m^2}.
\]

The total area, or volume per unit depth, decreases by

\[
1-\frac{A_f}{A_0}
=1-\frac{2475}{2700}
=0.08333,
\]

or approximately

\[
\boxed{8.33\%}.
\]

If this contraction produced a spatially uniform flat-surface drop, the equivalent final metal height would be

\[
h_f
=(5\times10^{-3})\frac{2475}{2700}
=4.5833\times10^{-3}\ \mathrm{m},
\]

corresponding to an average drop of

\[
5.0000-4.5833=0.4167\ \mathrm{mm}.
\]

The simulated free surface is not flat: contraction is concentrated into a central pipe while the solidified side regions remain higher. Consequently, \(0.4167\ \mathrm{mm}\) is an equivalent global-height change, not the predicted depth of the pipe at its center.

---

### 13. Hypothetical expansion calculation: Figure 13

The authors repeat the non-divergence-free calculation after reversing the densities:

\[
\rho^L=2700\ \mathrm{kg/m^3},
\qquad
\rho^S=2475\ \mathrm{kg/m^3}.
\]

Now the solid is less dense than the liquid. The metal must occupy a larger volume after solidification, and the low-Mach source gives

\[
\nabla\cdot\mathbf{u}>0
\]

inside the mushy region.

Figure 13 shows:

- solidification beginning from the two side walls;
- upward deformation of the gas-metal interface as the metal expands;
- a central protrusion developing while the remaining liquid channel closes; and
- a permanent raised feature after complete solidification.

The paper reports that this case also finishes solidifying at approximately \(t=0.25\ \mathrm{s}\).

A global mass-conservation estimate gives

\[
\frac{V_f}{V_0}
=\frac{\rho^L}{\rho^S}
=\frac{2700}{2475}
\approx1.09091.
\]

Thus, the hypothetical material expands by approximately

\[
\boxed{9.09\%}
\]

relative to its initial liquid volume. If the rise were uniform, the equivalent final height would be

\[
5\ \mathrm{mm}\times\frac{2700}{2475}
\approx5.455\ \mathrm{mm}.
\]

As in the shrinkage case, the actual deformation is localized, producing a protrusion rather than a uniformly raised flat surface.

---

### 14. Comparison of the three Section 6.3 calculations

| Calculation | Density ordering | Constraint during phase change | Sign of \(\nabla\cdot\mathbf{u}\) | Final free-surface behavior |
|---|---:|---:|---:|---|
| Artificial divergence-free case | \(\rho^S>\rho^L\) | Forced \(\nabla\cdot\mathbf{u}=0\) | Zero | Surface remains flat; no volume change |
| Physical shrinkage case | \(\rho^S>\rho^L\) | Low-Mach constraint | Negative | Surface caves inward; pipe defect forms |
| Hypothetical expansion case | \(\rho^S<\rho^L\) | Low-Mach constraint | Positive | Surface moves outward; protrusion forms |

This comparison shows that the sign and magnitude of phase-change-induced volume deformation are controlled by both the density difference and the rate of change of liquid fraction.

---

### 15. PCM mass conservation and grid refinement

Appendix C evaluates the PCM mass using

\[
m(t)
=
\int_\Omega
\left[
\rho^L(H\varphi)
+
\rho^S(H-H\varphi)
\right]
\,d\Omega.
\]

The relative mass error is

\[
\mathcal{E}(t)
=
\frac{|m(t)-m_0|}{m_0}\times100.
\]

For the shrinkage-pipe calculation, Figure C.16 shows that \(\mathcal{E}\) decreases as the grid is refined. The grids and time steps used in that study are:

| Grid | Time step |
|---:|---:|
| \(32\times32\) | \(10^{-4}\ \mathrm{s}\) |
| \(64\times64\) | \(10^{-4}\ \mathrm{s}\) |
| \(128\times128\) | \(10^{-5}\ \mathrm{s}\) |
| \(256\times256\) | \(10^{-5}\ \mathrm{s}\) |
| \(512\times512\) | \(10^{-6}\ \mathrm{s}\) |

The liquidus-solidus interval is

\[
\Delta T=T_{\mathrm{liq}}-T_{\mathrm{sol}}=10\ \mathrm{K}
\]

for all grids. At the finest resolutions, the paper reports an approximately

\[
\boxed{0.39\%}
\]

change in PCM mass. The error also decreases when the time step is reduced, although the corresponding data are not shown separately.

The reported error is not at machine precision. As discussed for the melting case, likely contributors within the authors' formulation are the non-conservative level-set representation of the gas-PCM interface and numerical evaluation of the low-Mach divergence source. The important trend is that the mass error decreases with refinement.

---

### 16. What Section 6.3 validates

Section 6.3 demonstrates that the method can:

- initiate solidification from thermally cooled boundaries;
- create and remove a liquid region until the metal becomes fully solid;
- couple solidification to a moving gas-metal free surface;
- produce contraction when \(\rho^S>\rho^L\);
- produce expansion when \(\rho^S<\rho^L\);
- form a shrinkage-pipe defect without prescribing its geometry;
- form a protrusion when the density ordering is reversed;
- switch off the volume source after phase change is complete;
- remain stationary after complete solidification; and
- reduce PCM-mass error under grid and time-step refinement.

Most importantly, it demonstrates a formulation-level requirement:

\[
\boxed{
\text{A variable-density phase-change model must constrain velocity consistently with its equation of state.}
}
\]

Different phase densities in the momentum and enthalpy equations are not enough if velocity is still forced to be divergence-free throughout the phase-change region.

---

### 17. What Section 6.3 does not establish

The example does not provide:

- comparison with a measured casting-pipe geometry;
- an analytical transient solution for the evolving free surface;
- a quantitative experimental validation of solidification time;
- a detailed study of contact angle or wetting at the gas-metal-wall junction;
- a parameter study of surface tension, gravity, wall cooling, or mold geometry; or
- a three-dimensional casting calculation.

It is therefore best viewed as a **mechanistic numerical benchmark** demonstrating the kinematic necessity of phase-change-induced velocity divergence, rather than as a complete predictive model of an industrial casting process.

---

### 18. Main conclusion of Section 6.3

The divergence-free calculation solidifies the metal without changing its external volume and therefore misses the casting defect. The low-Mach calculation correctly links the density change to the velocity field:

\[
\rho^S>\rho^L,
\quad
\frac{D\varphi}{Dt}<0
\quad\Longrightarrow\quad
\nabla\cdot\mathbf{u}<0,
\]

which causes the free surface to cave inward and form a shrinkage pipe.

Reversing the density ordering reverses the sign of the divergence and creates a protrusion instead. These paired tests show that the low-Mach enthalpy method captures both the direction and the macroscopic free-surface consequence of phase-change-induced volume change.

---

### Reference for Sections 6.2 and 6.3

R. Thirumalaisamy and A. P. S. Bhalla, "A low Mach enthalpy method to model non-isothermal gas-liquid-solid flows with melting and solidification," *International Journal of Multiphase Flow*, vol. 169, 104605, 2023. Sections 6.2-6.3, Figures 10-13, Table 1, and Appendix C.

