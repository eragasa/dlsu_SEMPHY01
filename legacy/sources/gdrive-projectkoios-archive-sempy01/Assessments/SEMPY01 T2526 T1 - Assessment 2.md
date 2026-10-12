## Long-Form Exercise 2 (Ch. 5–6 Focus)

**Submission:**  
- Submit a single PDF report (exported from Markdown/Word/LaTeX).  
- Submit your **code as a .zip file** (Python/Matlab/etc.) separately through Canvas.  
- All plots must be generated computationally. No hand-drawn plots.  
- Axes must be annotated with units.  

**Structure:**  
- Each section begins on a new page.  
- After each computational result, include a **Reflection/Discrepancy Analysis** subsection.  
- Final day: **Presentation** (10–15 min per student).  

---

## Task 1 — Drift and Diffusion Currents
- Implement the drift and diffusion current densities:
  $$
  J_{\text{drift}} = q n \mu E, \qquad 
  J_{\text{diff}} = q D \frac{dn}{dx}.
  $$
- Plot $J_{\text{drift}}$ vs. $E$ for electron and hole mobilities.  
- Plot $J_{\text{diff}}$ for a linear carrier gradient.  
- Verify the **Einstein relation**:
  $$
  \frac{D}{\mu} = \frac{k_BT}{q}.
  $$

---

## Task 2 — Mobility and Temperature Dependence
- Implement a model $\mu(T) \propto T^{-m}$, where $m$ reflects dominant scattering.  
- Compute and plot $\mu(T)$ over $200 \,\mathrm{K} \leq T \leq 600 \,\mathrm{K}$.  
- Compare electrons vs. holes in your assigned semiconductor.  
- Reflect on how scattering mechanisms (phonons, impurities) influence the result.  

---

## Task 3 — Nonequilibrium Excess Carriers
- Solve the continuity equation for excess carriers under steady-state illumination:  
  $$
  \frac{\partial \Delta n}{\partial t} = G - \frac{\Delta n}{\tau}.
  $$
- Plot $\Delta n(t)$ showing recombination with carrier lifetime $\tau$.  
- Explore at least two cases: **short lifetime** vs. **long lifetime**.  
- Compare your results with literature values of lifetimes.  

---

## Task 4 — Quasi-Fermi Levels
- For a forward-biased junction, compute quasi-Fermi levels $E_{Fn}$ and $E_{Fp}$.  
- Plot $E_{Fn}$ and $E_{Fp}$ relative to $E_c$ and $E_v$.  
- Discuss the physical meaning:  
  - Separation $\Delta E_F = qV_a$ corresponds to applied bias.  

---

## Task 5 — Experimental Comparison
- Find at least one dataset (mobility, diffusion constant, or lifetime) for your semiconductor.  
- Overlay with your computational results.  
- Identify one **agreement** and one **discrepancy**.  
- Suggest at least one reason for mismatch.  

---

## Reflection and Discrepancy Analysis (15%)
Each section must include a short reflection.  

### Template
1. **What did you expect?**  
   (Describe the anticipated trend or limit.)  

2. **What did you find?**  
   (Summarize the computed result.)  

3. **If there’s a mismatch, why?**  
   (Give at least one explanation in 2–3 sentences.)  

---

## Grading Breakdown
- Computational tasks and plots: **45%**  
- Literature comparison: **15%**  
- Reflection and discrepancy: **15%**  
- Presentation: **25%**  

---

## Report Format
Your report must be a single PDF with the following structure:

### 1. Introduction
- Identify your assigned semiconductor.  
- State why transport and recombination are important.  
- Outline the scope: drift, diffusion, mobility, lifetimes, quasi-Fermi levels.  

### 2. Methods
- Present equations used (with sources).  
- Note assumptions (e.g., constant effective mass, uniform doping).  
- Include short code snippets where relevant.  

### 3. Results
- Drift/diffusion current plots  
- $\mu(T)$ vs. $T$  
- $\Delta n(t)$ for two lifetimes  
- Quasi-Fermi level separation plot  
- Experimental comparison overlay  

### 4. Reflection and Discrepancy Analysis
- Follow the template for each result.  

### 5. Discussion
- Summarize insights: role of mobility, lifetimes, quasi-Fermi levels.  
- Compare results to literature.  
- Note limitations of your approach.  

### 6. Conclusion
- Restate your key findings.  
- Identify one thing you learned from coding.  
- State one open question about your material.  

### 7. References
- IEEE style.  
- Include at least **two sources** (textbook + journal).  
- Submit PDFs of cited references in your archive.  

---

## Midterm Presentation Format (10–15 min)
1. **Introduction** — Semiconductor and why it matters.  
2. **Methods** — Equations, assumptions, code structure.  
3. **Key Results** — Drift/diffusion, mobility, lifetime plots, quasi-Fermi levels.  
4. **Experimental Comparison** — Overlay + discussion.  
5. **Reflection** — One agreement, one discrepancy, one limitation.  
6. **Conclusion** — One next step.  
7. **Live Demo** — Re-run one plot (e.g., $\mu(T)$).  

---

## Academic Integrity
- References must be real and included as PDFs.  
- AI-generated or fabricated citations will be penalized.  
- You may use coding help (AI, StackOverflow, etc.) but must **explain your methods** in your own words.  
