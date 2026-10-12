[[SzeLee3Ed.pdf]] Section 2.05

#### Definition (Electron affinity)
Electron affinity $\chi$ is the difference between the conduction band minimum $E_c$ and the vacuum level $E_\text{vac}$.$$
q \chi= E_\text{vac}-E_c
$$
#### Definition (Work Function)
The work function $\psi_s$ is the energyt difference between the Fermi level $E_f$ and the vacuum level
$$
  q\phi_s=E_\text{vac}-E_f
$$
So we combine these two definitions
$$
  E_c-E_f=qV_s
$$
$$
q\phi_s=q(\chi+V_s)
$$
The barrier for electron emission is $q\chi$, not the work function.  Only the electrons with
$$
  E \ge E_\text{vac}=E_c+q\chi
$$
can escape.

![[Pasted image 20251118130435.png]]
- vacuum level is the highest reference point
- $E_c$ sits below by $q\chi$
- $E_f$ sits below the $E_c$ by $qV_s$

**Thermionic emmission requires the electrons in the conduction band tail whose energy exceeds the vacuum level by thermal excitation.**

If we are dealing with thermal excitation, then we have to get our information from statistical mechanics.

![[Pasted image 20251118130800.png]]
the maxwell boltzman distribution tail in the conduction band.

#### Count the emittable electrons

$n(E)$: density of states
$f(E)$: fermi occupation level

for non degenerate semiconductors
- $f(E) \approx \exp\!\left[-\frac{E - E_f}{k_B T}\right]$
- $n(E)=N_c \exp \left[ - \frac{E-E_c}{k_B t}\right]$
	- conduction band effective DOS

Substitute the lower limit $E_c + q\chi$:
$$
  n_\text{th}= N_c \exp\left[
    - \frac{E_c+q\chi-E_f}{k_B T}
  \right]
$$
Since $E_c - E_f = q V_s$
$$
  n_\text{th}= N_c \exp\left[
    - \frac{q(\chi+V_s)}{k_B T}
  \right]
$$
#### Themionic Emission Current Density
The current comes from the flux of electrons with sufficient velocity directed toward the surface

In the region, where we can use Maxwell Boltzman distribution
$$
  J_{\text{TE}}=A^*T^2\exp\left[
    -\frac{q \phi_s}{k_B T}
  \right]
$$
$A^*$: is the richardard constant which dependent of effective mass.
$$
A^*=\frac{4\pi q m^*k^2}{h^3}
$$
$\psi_s = \chi + V_s$ is the semiconductor work function