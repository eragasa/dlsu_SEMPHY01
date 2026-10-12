## The Intrinsic Semiconductor
- [[Neaman3Ed.pdf|Neaman, 3rd Ed, Section 3.1 - Charge Carrriers in Semiconductors]]

Semiconductor has 
- the bottom of the conduction band, $E_c$
- the top of the valence band, $E_v$
- valence band ($E<E_F$)
- conduction band ($E > E_C$)
- Fermi Energy, $E_v < E < E_c$
- At $T=0 \text{ K}$, valence fully occupied and the conduction band is empty
- As $T$ increases, more electron leave the the valence band and occupy the conduction band.
- The reason this happens can be modelled by the Fermi-Dirac distribution

Current the rate which charge flows
- Semiconductors have two types of charges we care about
	- "electrons" - electrons which are in the conduction band ($n$)
	- "holes" - "vacancies" left by electrons that are in the conduction band
- There are density of states in both the valence band and the conduction band
	- so a "hole" is essentially an unoccupied density of state in the valence band
	- an "electron" is a occupied density of state in the conduction band

### Equilibrium Distribution of Holes and Electrons (Intrinsic Conductor)

#### Distribution of electrons in the conduction band
- $n(E)$ - density of electron in conduction band
- $g_c(E)$ - density of states for the conduction band
- $f_F(E)$ - the fermi dirac distribution
$$
	n(E)=g_c(E)f_F(E)
$$
The total electron concentration per unit volume in the conduction band is found be integrating over the entire conduction band energy
$$
	n_0 = \int_{E_c}^\infty g_c(E)f_F(E)\mathrm{d}E
$$
#### Distribution of holes in the valence band
Since all probability must sum to 1.   If an electron is not in conduction band, it is the valence band.
- $p(E)$ - density of a hole in the valence band
- $g_v(E)$ - density of states for the conduction band
- $f_F(E)$ - the fermi dirac distribution
- $(1-f_F(E))$ is the probability that a state is not occupied by an electron.
$$
	p(E)=g_v(E)[1-f_F(E)]
$$
Likewise, the total hole concentration per unit volume is found by integration this function over the entire balance band energy
$$
	p_0 = \int_{-\infty}^{E_v}g_v(E)[1-f_F(E)]\mathrm{d}E
$$
![[Pasted image 20250916131236.png]]


### Effective mass of an eletron
The **effective mass** $m^*$ describes how electrons and holes respond to forces inside a crystal lattice.  It accounts for the fact that carriers move in a periodic potential, not free space.

The effective mass of an electron $m_n^*$ and effective mass of a hole $m_p^*$.  
### Band Structure Approximation
Recall that the dispersion relationship related $E(\mathbf{k})$ vs $\mathbf{k}$.  Near the band extrema, the dispersion can be approximated as parabolic:
$$
E(\mathbf{k}) \approx E_0 + \frac{\hbar^2 |\mathbf{k} - \mathbf{k}_0|^2}{2 m^\star}
$$

- $E_0$: band edge energy ($E_c$ or $E_v$)  
- $\mathbf{k}_0$: extremum point in reciprocal space  
- $m^\star$: effective mass  
##### Definition: Effective Mass
The effective mass is related to the curvature of the band:
$$
\frac{1}{m^\star} = \frac{1}{\hbar^2}\,\frac{\mathrm{d}^2 E}{\mathrm{d}k^2}
$$

### Thermal Equilibrium for Intrinsic Semiconductor
#### T.E. for electron in the conduction band for intrinsic semiconductor
If we do the integration of 
$$
	n_0 = \int_{E_c}^\infty g_c(E)f_F(E)\mathrm{d}E
$$
We can get the following result
$$
n_0 =2\left(\frac{2\pi m_n^* kT}{h^2}\right)\exp\left[\frac{-(E_c-E_F)}{kT}\right]
$$
##### Effective Density of States
$$
	N_c = 2\left(\frac{2\pi m_n^* kT}{h^2}\right)
$$
##### Electron concentration in conduction band
$$
	n_0=N_c\exp\left[\frac{-(E_c-E_F)}{kT}\right]
$$
#### T.E. for holes in the valence band for intrinsic semiconductor
Evaluating the the integral
$$
	p_0 = \int_{-\infty}^{E_v}g_v(E)[1-f_F(E)]\mathrm{d}E
$$
We get
$$
	p_0=2\left(\frac{2\pi m^*_p kT}{h^2} \right)\exp\left(\frac{-(E_F-E_v)}{kT}\right)
$$

##### Effective density of states in the valence band
$$
	N_v = 2\left(\frac{2\pi m^*_p kT}{h^2} \right)
$$
##### Hole concentration
$$
	p_0=N_v\exp\left(\frac{-(E_F-E_v)}{kT}\right)
$$
##### Example: Silicon
Calculate the probability that a state in the conduction band is occupied by an electron and calculate the thermal equilibrium electron concentration in silicon at $T=100K$
1. Assume the fermi energy is 0.25 eV below the conduction band.   The value of $N_c$ for silicon at $T=100\text{ K}$ is $N_c=2.8\times10^{19}\text{ }\mathrm{cm^{-3}}$.
2. Holes too!
## Intrinsic Carrier Concentration

- $(E_{F})_i$ - The Fermi energy level for the intrinsic semiconductor
- $n_i$ = electron concentration in the conduction band for the intrinsic semiconductor ("intrinsic electron concentration")
- $p_i$ = hole concentration in the valence band for the intrinsic semiconductor ("intrinsic hole concentration)

Since in the intrinsic semiconductor all electrons in the conduction band come from the valence band, then concentrations of holes and electrons must be the same
$$n_i=p_i$$
So we can express the equilibrium condition as
$$
  np=n^2
$$
$$
n_i^2 = N_c N_p \exp \left[ \frac{-(E_c - E_v)}{kT} \right]
$$

If we define the bandgap as 
$$
	E_g=E_c-E_v
$$
$$
\tag{Neaman 4.42}
n_i^2 = N_c N_p \exp \left[ \frac{-(E_g)}{kT} \right]
$$
- $N_c (T, m_c^*)$
- $N_p(T,m_v^*)$

##### Problem: 
Calculate the thermal equilibrium electron concentration in silicon for the case when $E_F=E$ and $T=300K$.

$\text{Ans: } 1.9 \times 10^19 \text{ }\mathrm{cm^{-3}}$
#### The Intrinsic Fermi-Level Position
For semiconductors, we kind always said that $E_{\text{midgap}}=\frac{1}{2}(E_c+E_v)$.   

Specifically calculate the fermi energy for the intrinsic semiconductor $(E_F)_i$ 
$$
	n_i=p_i
$$

$$
(E_{F})_i = E_{\text{midgap}} + \frac{3}{4} kT \ln \!\left( \frac{m_p^\star}{m_n^\star} \right)
$$
