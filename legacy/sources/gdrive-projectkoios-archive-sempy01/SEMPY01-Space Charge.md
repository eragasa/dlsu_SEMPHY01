[[SzeLee3Ed.pdf#page=83|Sze3rdEd, Section 2.7 - Space Charge Effect]]

# Motivation
In models we've studied, we assume charge neutrality
$$
  n+N_A^-\approx p+N_D^+
$$
so the internal electric field is small and carrier transport is dominated by drift-diffusion processes.  But this assumption breaks down in any device that requires carrier injection rather than pure doping.

Whenever push a contact to supply carriers into a semiconductor - LEDS, photodiotics, thin-film resistors, etc....
$$\begin{gather}
  n \ll N_D^+ \\
  p \gg N_A^-
  
\end{gather}$$
this makes the semiconductor temporarily non neutral
- injected charge creates it's own electric field (def. electric field)
- this field now controls how the carriers move.
- current-voltage relationship becomes non-linear

This injected charge is **space charge**.

## Space charge

#### Definition: Space charge
Space charge measures the _local imbalance_ between free carriers and ionized dopants. the space charge density is 
$$
  \rho =q(p-n+N_D^+N_A^-)
$$
where 
- $q$ — elementary charge
 - $n$ — electron concentration
 - $p$ — hole concentration
- $N_D^+$ — ionized donor concentration
- $N_A^-$ — ionized acceptor concentration
$\blacksquare$

#### Example: Neutral Region
- $n=N_D^+$
- $p =N_A^-$
- which means $\rho=0$
#### Example: $n$ type conductor
- $n \gg N_D^+$
- $p \approx N_A^- \approx 0$
- $\rho \approx -qn$
The injected carrier population becomes the soruce of electrostatic charge the determines the electric field through the Poisson equation.

## Drift Current Under Injection

For electron injection, the drift current
$$
  J =qnv
$$
Assume constant mobility
$$
  \nu=\mu E
$$
so
$$
J=qn\mu E
$$
Solving for $n$ relates $n$ to the local electric field
$$
\tag{1}
n=\frac{J}{q\mu E}
$$
## Poisson Equation with Injected Carriers
We have the assumption
- $n \gg N_D^+$
- $p \approx 0$
$$
\tag{2}
\frac{\mathrm{d}E}{\mathrm{d}x}
=\frac{\rho}{\epsilon_s}
=\frac{qn}{\epsilon_s}
$$
substitute Eqn (1) into Eqn (2)
$$
\tag{3}
\frac{\mathrm{d}E}{\mathrm{d}x}=\frac{J}{\epsilon_s\mu E}
$$
This is the differential euqation that governs SCLC (space charge limited current).   We are talking about SCLC because this is the simplest non-linear model for injected charges.

- If $J$ is fixed, the electric field must increase with distance $x$ to maintain the current
- **This leads to a nonuniform field profile inside the semiconductor**

## Integrate the Field Equation

Write Eqn (3) as 
$$
 E \mathrm{d}E=\frac{J}{\epsilon_s \mu}\mathrm{d}x
$$
we choose the boundary condition $E(0)=0$
$$
  E^2/2=\frac{J}{\epsilon_s \mu} x
$$

$$
E(x)=\sqrt{\frac{2J}{\epsilon_s \mu}}x^{1/2}
$$
This gives a square root field profile of the the SLDC

## Voltage-Current Relationship

$$
  V = \int_0^L E(x) \mathrm{d}x
$$
solving this integral
$$
  J=\frac{9}{8}\frac{\epsilon_s \mu}{L^2}V^2
$$
#### Definition: Mott-Gurney Law
space charge limited curent in the constant mobility regime scales as 
$$
  J \propto V^2
$$
## Physical Interprestation
- injected carriers product an internal space charge
- this space charge bends the bands (Sze Fig 21a)
- 
- ![[Pasted image 20251118140153.png]]
- the electric field becomes position dependent
- to sustain the curent, the voltage must rise faster linearly
- the mobility $\mu$ controls dirft,  higher $\mu$ means higher SCLC
- The thickness $L$ is critical, SCLC falls like $1/L^2$

This is why space charge is prominant in thin films, organic semiconductors, and III-V high-purity payers.

## Velocity Saturation Regime
Sze figure 3.21b
![[Pasted image 20251118140513.png]]
