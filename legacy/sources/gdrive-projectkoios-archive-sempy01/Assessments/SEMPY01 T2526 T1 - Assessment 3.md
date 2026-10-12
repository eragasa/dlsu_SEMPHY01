# Semiconductor Device-Level Modeling
#### SEMPY01 -Semiconductor Physics
#### AY 2025-2026 T1

**Format**  
- This is an **individual assessment**.  
- You must implement and analyze **all four canonical devices**:  
  1. PN junction diode  
  2. Bipolar junction transistor (BJT)  
  3. MOSFET  
  4. LED or photodiode  
- Device modeling must use material parameters drawn from Assessments 1 and 2.  
- Code must be your own, with clear plots and labeled axes.  
## Part A — PN Junction Diode
- Implement the diode equation:  
  $$
  I(V) = I_s\left(e^{qV/k_BT} - 1\right),
  $$
  where $I_s$ depends on $n_i$, diffusion coefficients, and diffusion lengths.  
- Plot the diode I–V curve over forward and reverse bias.  
- Compare behavior against a silicon diode baseline.  
## Part B — Bipolar Junction Transistor (BJT)
- Model a common-emitter configuration.  
- Show how emitter/base/collector doping asymmetry affects current gain $\beta$.  
- Plot $I_C$ vs $V_{CE}$ for multiple fixed $I_B$.  
- Discuss how your semiconductor’s $E_g$, $n_i$, and mobility affect transistor performance.  
## Part C — MOSFET
- Derive and implement drain current equations for the linear and saturation regions.  
- Plot $I_D$ vs $V_{DS}$ for at least three different $V_{GS}$ values.  
- Comment on threshold voltage and channel mobility in your chosen material.  
## Part D — LED or Photodiode
- **LED case**: use $E_g$ to estimate emission wavelength; calculate recombination rate and photon output.  
- **Photodiode case**: model photocurrent vs. photon energy; show spectral cutoff at $E_g$.  
- Include one literature band diagram or spectral curve for comparison.  
## Reflection and Discrepancy Analysis
Each section must include a short reflection:  
1. What did you expect?  
2. What did you find?  
3. If there was a mismatch, why?  
## Grading Breakdown
- PN junction diode: **25%**  
- BJT: **25%**  
- MOSFET: **25%**  
- LED/Photodiode: **15%**  
- Reflection & discrepancy analysis: **10%**  
