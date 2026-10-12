## Long-Form Midterm (Chapter 4 Focus)

**Submission:**  
- Submit a single PDF report (exported from Markdown/Word/LaTeX).  
- Submit your **code as a .zip file** (Python/Matlab/etc.) separately through Canvas.  
- All plots must be generated computationally (Python/Matplotlib, Matlab, etc.). No hand-drawn plots.  
- Annotate axes clearly with units.  
**Structure:**  
- Each section must begin on a new page.  
- After each plot or computational result, include a **Reflection/Discrepancy Analysis** subsection (see rubric below).  
- Final day: **Presentation** (10–15 min per student).  

## Task 1 — Intrinsic Carrier Concentration $n_i(T)$
- Derive the expression for $n_i(T)$ for your assigned semiconductor.  
- Implement the computation over a temperature range (200 K – 600 K).  
- Plot $n_i(T)$ on both linear and logarithmic scales.  
- Compare with literature values where possible.  

## Task 2 — Fermi Level $E_F(T)$
- Calculate the Fermi level as a function of temperature for your semiconductor.  
- Include intrinsic, doped-$n$, and doped-$p$ cases.  
- Plot $E_F(T)$ with respect to conduction and valence band edges.  
- Discuss limits ($T \to 0$, $T \to \infty$).  

## Task 3 — Incomplete Ionization
- Incorporate incomplete ionization of dopants into your model.  
- Show the freeze-out, extrinsic, and intrinsic regimes clearly.  
- Plot carrier concentration vs. $T$ including incomplete ionization.  
## Task 4 — Literature Review and Band Diagram
- Find at least **two references** (journal, textbook, or reliable database).  
- Identify one real **device application** for your semiconductor.  
- Explain why its properties (band gap, $n_i$, mobility, etc.) make it suitable.  
- Include a **band diagram** or **band structure figure** from the literature (with citation).  
## Task 5 — Experimental Comparison
- Locate at least one piece of **experimental data** for your semiconductor ($n_i(T)$, $E_g(T)$, mobility, etc.).  
- Overlay with your computational result where possible.  
- Highlight at least one **agreement** and one **discrepancy**, and explain.  
## Reflection and Discrepancy Analysis (15%)
Each section of the midterm must include a **Reflection/Discrepancy** subsection.  
This is graded separately from correctness.  

### Grading Breakdown
- **Method Explanation (5%)**  
  - Briefly describe how you set up your calculation or code (what formula, what loop/solver, what assumptions).  
  - Credit is for clarity, not correctness.  

- **Expected vs. Actual (5%)**  
  - State what you expected the trend/curve to look like (e.g., “$n_i(T)$ should increase exponentially with $T$”).  
  - Compare this expectation with your computed plot.  
  - Highlight any mismatches.  

- **Interpretation (5%)**  
  - Offer at least one possible reason for discrepancy.  
  - This can be physics-based (wrong assumption), math-based (units, constants), or computational (plotting/log scale).  

### Reflection Template
After each plot or result, answer:

1. **What did you expect?**  
   (Describe the anticipated behavior: shape, limits, trends.)  

2. **What did you find?**  
   (Summarize the actual computed result or plot.)  

3. **If there’s a mismatch, why?**  
   (Give at least one explanation in 2–3 sentences.)  
## Long-Form Midterm Report Format
Your report must be submitted as a single PDF.  
Each section should start on a **new page**.  
Figures must be labeled and referenced in the text.  
Keep the write-up concise and technical — avoid filler text.
## 1. Introduction
- Identify your assigned semiconductor.  
- State why it is important in semiconductor physics.  
- Summarize the **scope of your project** (carrier concentration, Fermi level, ionization, comparison with experiment).  
- Include references to textbooks or papers (IEEE style).  
## 2. Methods
- Clearly state the **equations used** (with sources).  
- Describe **assumptions** (intrinsic vs. extrinsic, complete vs. incomplete ionization).  
- Provide short **code snippets** if necessary to illustrate your approach (not the full code).  
- Mention the constants and parameters you used ($k_B$, $\hbar$, $E_g$, etc.).  
## 3. Results
- Present plots in order:
  1. $n_i(T)$ vs. $T$  
  2. $E_F(T)$ for intrinsic and doped cases  
  3. Carrier concentration with incomplete ionization  
  4. Literature band diagram (with citation)  
  5. Comparison with experimental data (overlay if possible)  

- Each figure should have a **caption** and be referenced in the text.  
## 4. Reflection and Discrepancy Analysis
For each result:
1. What did you expect?  
2. What did you find?  
3. If there was a mismatch, why?  
## 5. Discussion
- Summarize key physics insights (e.g., why $n_i$ rises with $T$, why incomplete ionization matters).  
- Compare with literature: how close were your results?  
- Comment on limitations of your model (constants, approximations, missing effects).  
## 6. Applications
- From your literature review:  
  - Describe one device that uses your semiconductor.  
  - Explain (briefly) how the material properties make it suitable.  
## 7. Conclusion
- Restate your main findings.  
- Identify one thing you learned from coding.  
- Identify one open question about your material.  
## References
- Use IEEE style citations.  
- At least **two journal/textbook references** are required.  

## Grading Breakdown
- Computational tasks and plots: **45%**  
- Literature review and band diagram: **15%**  
- Experimental comparison: **10%**  
- Reflection and discrepancy analysis: **15%**  
- Presentation: **15%**  

## Submission Requirements
1. **Final Report (PDF)**  
   - Structured according to the Report Format.  
   - Figures must be labeled, referenced in the text, and explained.  

You must submit a **.zip archive** to Canvas containing:
1. **Code Folder**  
	- All scripts used to generate your results.  
	- Code must run without modification and reproduce the plots.  
	- Organize clearly (e.g., `ni_vs_T.py`, `fermi_level_doping.py`, `incomplete_ionization.py`).
2. **References Folder**  
   - PDFs of **every reference you cite** in your report.  
   - References from textbooks (can just be the page number)
   - At least one peer-reviewed journal article or reliable source.  
   - Files must be named consistently (google scholar name from .bib). 
   - the relevant portions you used should be highlighted.
## Midterm Presentation Format
**Time Limit:** 10–15 minutes per student  
**Deliverables:** Slides (PDF or PPTX) + Live demonstration of code/plots  

### Required Sections
1. **Introduction (1–2 min)**
   - Assigned semiconductor + why it’s important
   - Quick overview of band gap, typical applications

2. **Methods (2–3 min)**
   - Show your code workflow (equations implemented, parameters used)
   - Explain one loop/function/solver in detail
   - Justify any assumptions (e.g., constant effective mass)

3. **Key Results (4–5 min)**
   - $n_i(T)$ (linear + log scale)  
   - Fermi level vs. $T$ (intrinsic, n-type, p-type)  
   - Incomplete ionization plot  
   - At least one comparison with **Si results** (sanity check)

4. **Literature Review & Band Diagram (2–3 min)**
   - Show one band diagram or band structure figure from a cited source  
   - Briefly summarize at least one application for your semiconductor  
   - Connect properties (band gap, mobility, etc.) to application

5. **Experimental Comparison & Reflection (2–3 min)**
   - Overlay one experimental dataset with your results  
   - Identify at least one agreement and one discrepancy  
   - Provide a physical/computational explanation

6. **Conclusion (1 min)**
   - Summarize your main finding(s)  
   - State one limitation of your model  
   - Suggest one next step (e.g., adding temperature dependence of $E_g$)

---

### Live Demonstration
- Open your code and re-run one plot (e.g., $n_i(T)$ for Si).  
- Be ready to explain *what part of the code controls what part of the physics*.  

### Q&A (2–3 min)
- Instructor and classmates may ask probing questions:
  - “Why does $n_i$ increase with $T$?”  
  - “Where does the band gap enter your code?”  
  - “What changes if you double the doping concentration?”  

## Academic Integrity
- A citation in your report without a corresponding PDF in your archive will not be counted.  
- AI-generated or fabricated references will result in penalties.  
- All references must be **actual documents you read**, not generic Google Scholar links.  
