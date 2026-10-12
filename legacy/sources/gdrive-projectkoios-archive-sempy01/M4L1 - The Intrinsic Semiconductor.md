Learning Objectives (Neaman, Ch 4)
- Derive the equations for thermal equilibrium of electrons and holes in terms of Fermi energy
- Derive the equations for intrinsic carrier concentration
- Derive the expression for the intrinsic Fermi level


# Conceptual Derivations — Intrinsic Semiconductors

## 1) Thermal equilibrium of electrons and holes in terms of the Fermi energy

**Prompt.**  
Explain why a semiconductor at thermal equilibrium is described by a **single, spatially constant Fermi level** $E_F$. Use this to relate the electron and hole concentrations to $E_F$.

**Assumptions.**  
- Equilibrium (no external drive): quasi-Fermi levels coincide, $E_{Fn}=E_{Fp}=E_F$.  
- Nondegenerate limit for simplicity (Boltzmann tail near band edges).  

**Key relations to invoke.**  
- Fermi–Dirac occupancy: $f(E)=\dfrac{1}{1+\exp[(E-E_F)/kT]}$.  
- Near edges (Boltzmann):  
  $$
  n \approx N_c \exp\!\left(-\frac{E_c-E_F}{kT}\right),\qquad
  p \approx N_v \exp\!\left(-\frac{E_F-E_v}{kT}\right).
  $$

**Targets.**  
- Argue **why** $E_F$ must be constant at equilibrium (no carrier or energy flow when $\nabla E_F=0$).  
- Express $n(E_F)$ and $p(E_F)$ as above and state the **law of mass action**: $np = n_i^2$ (independent of doping at equilibrium).

**Minimal outline.**  
Equilibrium $\Rightarrow$ one $E_F$; insert $E_F$ into FD occupancy; integrate DOS $\Rightarrow$ exponential forms for $n,p$; multiply to get $np$ constant.

## Intrinsic carrier concentration $n_i$

**Prompt.**  
Starting from the band-edge density-of-states (DOS) and Fermi–Dirac statistics, derive the intrinsic carrier concentration $n_i$.

**Assumptions.**  
- Intrinsic material: $n=p=n_i$.  
- Nondegenerate limit near edges.  
- Parabolic bands $\Rightarrow$ effective masses $m_n^\star, m_p^\star$.

**Key relations to invoke.**  
- Effective DOS:
$$
  N_c = 2\!\left(\frac{2\pi m_n^\star kT}{h^2}\right)^{\!\!3/2},\qquad
  N_v = 2\!\left(\frac{2\pi m_p^\star kT}{h^2}\right)^{\!\!3/2}.
  $$
- Boltzmann forms (from Q1):
$$
  n = N_c e^{-(E_c-E_F)/kT},\quad
  p = N_v e^{-(E_F-E_v)/kT}.
$$

**Targets.**  
- Show
  $$
  n_i^2 = N_c N_v \exp\!\left(-\frac{E_g}{kT}\right),\qquad E_g=E_c-E_v.
  $$
- Hence
  $$
  n_i = \sqrt{N_c N_v}\;\exp\!\left(-\frac{E_g}{2kT}\right).
  $$

**Minimal outline.**  
Use $n_i=p_i$ and $np=n_i^2$ with the Boltzmann forms; eliminate $E_F$ to obtain $n_i$ in terms of $N_c,N_v,E_g$.


---

## 3) Intrinsic Fermi level $E_{Fi}$ (position relative to midgap)

Derive the location of the **intrinsic Fermi level** $E_{Fi}$ relative to midgap, and show how it depends on the effective masses.

**Assumptions.**  
- Intrinsic: $n=p=n_i$.  
- Nondegenerate, parabolic bands with $m_n^\star, m_p^\star$.  
- Use $N_c,N_v$ from Q2.

**Key relations to invoke.**  
From Q1:
$$
n = N_c e^{-(E_c-E_{Fi})/kT},\quad
p = N_v e^{-(E_{Fi}-E_v)/kT}.
$$
Set $n=p$ (intrinsic).

**Targets.**  
- Solve $N_c e^{-(E_c-E_{Fi})/kT} = N_v e^{-(E_{Fi}-E_v)/kT}$ for $E_{Fi}$ to get:
  $$
  E_{Fi} = \frac{E_c+E_v}{2} + \frac{kT}{2}\ln\!\left(\frac{N_v}{N_c}\right).
  $$
- Substitute $N_c,N_v$ to obtain the **mass-dependent shift**:
  $$
  E_{Fi} - E_{\text{midgap}} 
  = \frac{3}{4}kT \ln\!\left(\frac{m_p^\star}{m_n^\star}\right),
  \qquad E_{\text{midgap}}=\frac{E_c+E_v}{2}.
  $$
- Interpret the sign: if $m_p^\star>m_n^\star$ then $E_{Fi}$ lies **above** midgap; if $m_p^\star<m_n^\star$ then **below**.

**Minimal outline.**  
Equate $n$ and $p$; isolate $E_{Fi}$; express $\ln(N_v/N_c)$ via $m^\star$; simplify to the $3/4\,kT$ result.
