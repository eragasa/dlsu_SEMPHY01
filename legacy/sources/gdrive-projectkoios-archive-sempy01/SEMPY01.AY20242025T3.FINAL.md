# SEMPY01 - FINAL ASSESSMENT
## Submission Requirements
- Submit a **single PDF report** (exported from Markdown/Word/LaTeX).  
- Submit your **code** as a `.zip` archive.  
- All plots must be generated computationally.  
- Axes must include **units**.  
- Each section starts on a **new page**.  
- Every result must include a **Reflection / Discrepancy Analysis** subsection.  
- Final deliverable includes a **10–15 minute presentation**.

# Task 1 — Drift and Diffusion Currents

Implement the current densities  
$$
J_{\mathrm{drift}} = q n \mu E, \qquad
J_{\mathrm{diff}} = q D \frac{dn}{dx}.
$$

Requirements:
1. Plot $J_{\mathrm{drift}}$ vs. $E$ for **electrons** and **holes**.  
2. Assume a linear carrier gradient and plot $J_{\mathrm{diff}}$.  
3. Verify the Einstein relation:
$$
\frac{D}{\mu} = \frac{k_B T}{q}.
$$

---

# Task 2 — Mobility and Temperature Dependence

Use the model  
$$
\mu(T) \propto T^{-m},
$$
where $m$ reflects the dominant scattering mechanism.

Requirements:
1. Compute and plot $\mu(T)$ for $200\,\mathrm{K} \le T \le 600\,\mathrm{K}$.  
2. Compare electron vs. hole mobility.  
3. Discuss phonon vs. impurity scattering.

---

# Task 3 — Nonequilibrium Excess Carriers

Solve the continuity equation  
$$
\frac{\partial \Delta n}{\partial t} = G - \frac{\Delta n}{\tau}.
$$

Requirements:
1. Plot $\Delta n(t)$ for at least two lifetimes (short vs. long).  
2. Compare with literature lifetimes for your semiconductor.

---

# Task 4 — Quasi-Fermi Levels

For a forward-biased $p$–$n$ junction, compute and plot the quasi-Fermi levels.

Requirements:
1. Plot $E_{Fn}$ and $E_{Fp}$ relative to $E_c$ and $E_v$.  
2. Explain the meaning of the separation:
$$
\Delta E_F = E_{Fn} - E_{Fp} = qV_a.
$$

---

# Task 5 — Experimental Comparison

Find one real dataset for your semiconductor (mobility $\mu$, diffusion constant $D$, or lifetime $\tau$).

Requirements:
1. Overlay the experimental dataset with your computed results.  
2. Identify one **agreement** and one **discrepancy**.  
3. Provide at least one **physical or computational** explanation.

---

# Task 6 — PN Junction Diode

Implement the diode equation  
$$
I(V) = I_s\left(e^{qV/k_B T} - 1\right).
$$

Requirements:
1. Plot $I$–$V$ for forward and reverse bias.  
2. Compare with a **baseline silicon diode**.  
3. Comment on effects of band gap, intrinsic concentration, and carrier properties.

# Report Structure (PDF)
## 1. Introduction
- Identify your semiconductor.  
- Explain why transport, mobility, recombination, and quasi-Fermi levels matter.  
- Outline Tasks 1–6.
## 2. Methods
- Present equations used.  
- State assumptions (uniform doping, constant effective mass, low-level injection, etc.).  
- Include short code snippets where relevant.
## 3. Results
Include:
- Drift current plot  
- Diffusion current plot  
- Einstein relation check  
- $\mu(T)$ plot  
- $\Delta n(t)$ for two lifetimes  
- Quasi-Fermi level plot  
- Experimental overlay  
- Diode I–V curve  
## 4. Discussion
- Summarize insights about mobility, lifetime, transport, quasi-Fermi levels, diode behavior.  
- Compare with literature.  
- State limitations.
## 5. Conclusion
- Restate key findings.  
- State one thing learned from coding.  
- Pose one open question about your semiconductor.
## 6. References
- IEEE style.  
- At least two references (textbook + journal).  
- Include PDFs of all references in your submission archive.
---
# Presentation Format (10–15 minutes)

1. **Introduction** — Semiconductor + motivation  
2. **Methods** — Equations, assumptions, code structure  (A FLOW CHART OF COMPUTATION IS FINE.  no need to do a explanation of methods)
3. **Results** — Drift/diffusion, mobility, lifetime, quasi-Fermi levels, diode I–V  
4. **Experimental Comparison** — Overlay + commentary  
5. **Reflection** — Agreement, discrepancy, limitation  
6. **Conclusion** — One next step  

