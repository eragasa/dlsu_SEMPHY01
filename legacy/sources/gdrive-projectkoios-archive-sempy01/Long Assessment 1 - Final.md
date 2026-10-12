# Long-Form Midterm Research Report  
**Course:** Semiconductor Physics  
**Focus:** Carrier Statistics and Energy Band Modeling
**Format:** Research Report + Code Demonstration

---

## Overview

This midterm is not a traditional problem set.  
It is a research-style investigation into carrier statistics, Fermi level behavior, and temperature-dependent properties of semiconductors.  

You will computationally and analytically model:
1. The **intrinsic carrier concentration** $n_i(T)$  
2. The **temperature-dependent Fermi level** $E_F(T)$  
3. **Incomplete ionization** and **freeze-out behavior**  
4. **Band structure interpretation** and **experimental validation**

Your final submission will combine analytical derivations, computational results, literature comparison, and scientific reflection.  Treat this as your first **mini research paper** in semiconductor physics.

---

## Deliverables

### (1) **Final Report (PDF)**
Exported from Markdown, Word, or LaTeX — must include:
- Equations, plots, and code explanations  
- Clearly labeled and referenced figures  
- Structured as a research article (see below)

### (2) **Code Archive (.zip)**
Include all scripts (Python/Matlab/etc.) used to generate figures:
- `ni_vs_T.py` – computes and plots $n_i(T)$  
- `fermi_level_T.py` – computes $E_F(T)$ for intrinsic, n-type, and p-type  
- `incomplete_ionization.py` – models $N_D^+$, $N_A^-$ vs. $T$  
- `compare_experiment.py` – overlays experimental and computational data  

Your code must reproduce all plots in your report when re-run.

### (3) **Reference Folder**
Include PDFs of every reference you cite (journal or textbook).  
Each file should be named by its BibTeX key and have the relevant section **highlighted**.  

(4) 5-10 minute presentation
 ---

## Report Structure

### 0. Abstract

### 1. Introduction

- Identify your assigned semiconductor (e.g., Si, Ge, GaAs, InP).  
- Explain its importance: typical band gap, applications, and **doping characteristics**.  
- Summarize the **objectives** of your study:  
  - Modeling $n_i(T)$  
  - Determining $E_F(T)$ for intrinsic and doped cases  
  - Exploring incomplete ionization and compensation  
  - Comparing with experimental or literature data

### **2. Methods**
This section replaces the “solution steps” of a problem set.  
Write as a formal methods section in a research paper.

Include:
1. **Governing Equations**
   - Intrinsic carrier concentration (Sze2014 Eq 21.14):
$$
 \tag{1}
 n_i(T) = \sqrt{N_C(T)\,N_V(T)}\,
 \exp\!\left(-\frac{E_g(T)}{2kT}\right)
$$
     with
 $$
 E_g(T) = E_g(0) - \frac{\alpha T^2}{T+\beta}.
 $$
   - Effective density of states:
 $$
 N_C(T) = 2\left(\frac{2\pi m_e^* kT}{h^2}\right)^{3/2}, \qquad
 N_V(T) = 2\left(\frac{2\pi m_h^* kT}{h^2}\right)^{3/2}.
 $$
   - Fermi level position:
 $$
 E_F - E_i = kT\ln\!\left(\frac{n_0}{n_i}\right)
 $$

- Incomplete ionization of donors:
 $$
 N_D^+ = \frac{N_D}{1+\exp\!\left[\frac{E_F-E_D}{kT}\right]}.
 $$

- Charge neutrality:
 $$
 n + N_A^- = p + N_D^+.
 $$

2. **Assumptions**
   - Non-degenerate statistics unless $E_F$ approaches $E_C$.
   - Full ionization at high temperature.
   - Partial ionization at cryogenic temperatures.
   - Effective masses are temperature-independent.

3. **Computational Implementation**
   - Mention constants ($k_B$, $q$, $\hbar$) and physical parameters.
   - Describe your workflow:
     - Temperature array: 200–600 K  
     - For each $T$, compute $E_g(T)$, $n_i(T)$, $E_F(T)$.
     - Use iterative or root-solving methods where needed.
   - Include a small code snippet (or psuedocode/not full script) to illustrate your implementation.
``` python
for i in some_array:
  k = get_Nd(i) # implements Eq(1)
```

---

### 3. Results

#### **3.1 Intrinsic Carrier Concentration $n_i(T)$**

- Compute and plot $n_i(T)$ (200–600 K).
- Show both **linear** and **logarithmic** scales.
- Overlay or compare with **literature data** (Sze, Neamen, etc.).

#### **3.2 Fermi Level $E_F(T)$**

- Compute $E_F(T)$ for intrinsic, $n$-type, and $p$-type cases. ( should probably put them on the same graph)
- Plot $E_F(T)$ relative to $E_C$ and $E_V$:
  $$
  E_C - E_i = \frac{E_g(T)}{2} - \frac{3}{4}kT\ln\!\left(\frac{m_h^*}{m_e^*}\right).
  $$
- Discuss low- and high-temperature limits:
  - $T \to 0$ (freeze-out)
  - $T \to \infty$ (intrinsic behavior)

#### **3.3 Incomplete Ionization and Freeze-Out**

- Incorporate $N_D^+$ and $N_A^-$ using Fermi–Dirac occupation.
- Plot **carrier concentration vs. temperature** for 0–600 K.
- Identify **three regimes**:
  1. Freeze-out (low $T$)
  2. Extrinsic (midrange)
  3. Intrinsic (high $T$)

#### **3.4 Literature and Band Diagram**

- Include at least **one band diagram or band structure** from a cited source.
- Use $E_i = 0$ as energy reference and label $E_C$, $E_V$, $E_F$.
- Summarize why your semiconductor’s **band gap** and **mobility** make it useful for a real device.

#### **3.5 Experimental Comparison**

- Find one dataset for $n_i(T)$, $E_g(T)$, or mobility from a journal or database (e.g., NIST, Ioffe, SemiconData).  
- Overlay your computational curve and comment:
  - Where do they agree?  
  - Where do they diverge?  
  - Possible causes (bandgap model, effective masses, degeneracy).

---

### **4. Reflection and Discrepancy Analysis**

**Reflection Block** directly after the figure or result.  
Use this structured template:

**Reflection Template**

1. **What did you expect?**  
   Describe your prior expectation or analytical prediction.  

2. **What did you find?**  
   Summarize the observed or computed trend.  

3. **If there’s a mismatch, why?**  
   Provide at least one plausible reason (physical, mathematical, or computational).

---

### **5. Discussion**

- Summarize physical interpretation:
  - Why does $n_i$ rise exponentially with $T$?
  - Why does $E_F$ shift toward $E_i$ at high $T$?
  - How does incomplete ionization manifest at low $T$?
- Comment on **degeneracy** and when Boltzmann approximation fails.
- Evaluate consistency with literature and experimental results.
- Identify **limitations** (e.g., neglect of impurity band formation at $N_D>10^{19}$ cm⁻³).

---

### **6. Applications**

- Identify a **device or technology** that uses your semiconductor (e.g., Si in CMOS, GaAs in LEDs).  
- Explain how specific properties — band gap, $n_i$, mobility, or incomplete ionization — determine its suitability.

---

### **7. Conclusion**

- Summarize major findings numerically:
  - Example: “For Si, $E_F - E_i = 0.41$ eV at $N_D = 10^{17}\,\text{cm}^{-3}$ and 300 K.”
- Identify one key concept learned through coding or model verification.
- Propose one follow-up question or possible model improvement.

---

### **8. References**

Use **IEEE citation style** throughout.  
At least **two peer-reviewed or textbook references** are required, such as:
- S. M. Sze and K. K. Ng, *Physics of Semiconductor Devices*, 3rd ed.  
- D. A. Neamen, *Semiconductor Physics and Devices*, 4th ed.  
- Journal sources (e.g., *J. Appl. Phys.*, *IEEE Trans. Electron Devices*).

Include corresponding PDF copies (highlighted sections) in your reference folder.

---

## Grading Rubric

| Component                             | Description                                                               | Weight |
| ------------------------------------- | ------------------------------------------------------------------------- | ------ |
| **Computational Modeling**            | Accuracy and reproducibility of $n_i(T)$, $E_F(T)$, incomplete ionization | 45%    |
| **Literature Review + Band Diagram**  | Proper citation, clarity, relevance                                       | 15%    |
| **Experimental Comparison**           | Data overlay, agreement/discrepancy discussion                            | 10%    |
| **Reflection & Discrepancy Analysis** | Conceptual insight, diagnostic reasoning                                  | 15%    |
| **Presentation & Code Demo**          | Clarity, reproducibility, command of material                             | 15%    |

---

## Midterm Presentation (10–15 Minutes)

Each student must present slides (PDF or PPTX) and demonstrate one code plot live.

**Required Sections:**
1. **Introduction:** Semiconductor overview and motivation  
2. **Methods:** Equations, model, assumptions, and computational flow  
3. **Results:** Plots for $n_i(T)$, $E_F(T)$, and incomplete ionization  
4. **Literature Review:** Band diagram and one device application  
5. **Experimental Comparison:** Overlay and discussion  
6. **Conclusion:** Summary + future work

**Q&A Focus Areas:**
- “Why does $n_i$ increase with $T$?”  
- “What physical parameter controls $E_F$?”  
- “Where does incomplete ionization appear in your equations?”  
- “How would degeneracy affect your results?”

---

## Academic Integrity

- Every cited reference must be backed by an actual PDF in your submission archive.  
- Generative or fabricated references will receive zero credit.  
- All computational results must be your own.  
- Collaboration may occur in *discussion*, not in *code sharing*.

---

**End of Midterm Research Brief**


Bibliography:
(Sze2021) - Sze, Simon M., Yiming Li, and Kwok K. Ng. _Physics of semiconductor devices_. 3rd Ed. John wiley & sons, 2014.