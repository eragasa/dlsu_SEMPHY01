## The Extrinsic Semiconductor

**Extrinsic semiconductor:**. A semi conductor that has a controlled amount of dopant atoms (either $n$-type donors or $p$-type acceptors)

![[Pasted image 20250919123326.png]]
Silicon: bonds with 4 neighbors
### $n$-type doping
Phosphorus: Group V elements: five valence electrons
- creates loosely bound electron
- donor impurityelectron
$E_C$: conduction band energy
$E_D$: donor level energy.   The discrete energy introduced by an $n$type dopant atom.   It lies just below the conduction band edge $E_C$.  When ionized, the donor contributes an electron to the electron band

![[Pasted image 20250919104915.png]]

### $p$-type doping

Boron: Group III element
$E_V$: valence band energy
$E_A$: acceptor level energy.  The discrete energy level introduced by a $p$-type dopant atom.  It lies just above the valence band $E_V$.  When ionized, the acceptor captures the electron, leaving behind a mobile hole in the valence vand



![[Pasted image 20250919105104.png]]

## Ionization Energy

# Ionization Energy in Semiconductors

Ionization energy in semiconductor physics usually refers to the **activation energy** required
to free an electron (or hole) bound to a dopant atom into the conduction (or valence) band.

---

## 1. Dopant Ionization Energy 
(Hydrogenic Model)

**Bohr Model**
The total energy of the orbiting electron is
$$ E = T + V $$
$$ T = \frac{1}{2}m^* v^2 $$
$$V = -\frac{e}{4\pi\epsilon r_n}=\frac{-m^*e^4}{(n\hbar)^2(4\pi\epsilon)^2}$$$$

A donor electron is bound to the positively charged dopant core in a manner similar to the
hydrogen atom, but modified by:

- The **effective mass** of the carrier $m^*$  
- The **dielectric constant** of the semiconductor $\varepsilon_r$

Effective Rydberg energy:
$$
E_R^* = \frac{m^*}{m_e} \cdot \frac{13.6 \,\text{eV}}{\varepsilon_r^2}
$$

Effective Bohr radius:
$$
a_B^* = \frac{\varepsilon_r}{(m^*/m_e)} \, a_B
\qquad (a_B = 0.529 \,\text{\AA})
$$


- **Donor level**: $E_C - E_D \approx E_R^*$ (uses electron mass $m_e^*$)  
- **Acceptor level**: $E_A - E_V \approx E_R^*$ (uses hole mass $m_h^*$)

---

### Examples

- **Silicon (Si):** $\varepsilon_r = 11.7$, $m_e^* \approx 0.26 m_e$  
  $E_D \approx 30$–$50$ meV

- **Germanium (Ge):** $\varepsilon_r = 16.0$, $m_e^* \approx 0.12 m_e$  
  $E_D \approx 10$ meV

- **Gallium arsenide (GaAs):** $\varepsilon_r = 12.9$, $m_e^* \approx 0.067 m_e$  
  $E_D \approx 5$–$6$ meV

---

## 2. Temperature Dependence (Freeze-Out → Extrinsic → Intrinsic)

At low $T$, dopants are **incompletely ionized** (freeze-out). As $T$ rises, all dopants ionize
(extrinsic regime), and at high $T$, intrinsic carriers dominate.

Fraction of ionized donors:

$$
\frac{N_D^+}{N_D} =
\frac{g_D \exp \!\left(\frac{E_D - E_F}{kT}\right)}{1 + g_D \exp \!\left(\frac{E_D - E_F}{kT}\right)}
$$

Fraction of ionized acceptors:

$$
\frac{N_A^-}{N_A} =
\frac{1}{1 + g_A \exp \!\left(\frac{E_F - E_A}{kT}\right)}
$$

where $g_D \approx 2$, $g_A \approx 4$ are degeneracy factors.

---

## 3. Extracting Ionization Energy from Experiment

In the freeze-out regime:

$$
n(T) \;\propto\; \exp\!\left(-\frac{E_C - E_D}{kT}\right)
$$

So an Arrhenius plot of $\ln n$ vs. $1/T$ has slope:

$$
\text{slope} = -\frac{E_C - E_D}{k}
$$

This gives the ionization energy.

---

## 4. Alternative Meaning: Pair Creation Energy

In detector physics, *ionization energy* can mean the **mean energy to create one electron–hole
pair**, which is several times the band gap:

- Si: ~3.6 eV per pair  
- Ge: ~2.9 eV  
- GaAs: ~4.2 eV

---

## Summary

- **Ionization energy (dopants):** meV scale, set by effective mass and dielectric constant.  
- **Incomplete ionization:** important at low $T$ (freeze-out).  
- **Experimental extraction:** Arrhenius analysis of $n(T)$ or $\sigma(T)$.  
- **Detector definition:** energy per $e$–$h$ pair, ~3–5 eV depending on material.
