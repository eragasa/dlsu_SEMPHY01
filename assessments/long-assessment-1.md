# Long Assessment 1: From Bloch Bands to Equilibrium Semiconductor Carrier Statistics

## Purpose and scope

In this assessment, you will construct and validate a computational model that connects microscopic band dispersion to the equilibrium electron and hole populations of a bulk semiconductor. You will begin with a numerical `Bloch2D` model, reduce its band edges to an effective-mass representation, and investigate when that reduction reproduces the local curvature and density of states of the parent calculation. You will then build a three-dimensional intrinsic-semiconductor model, apply it to a real material chosen from the literature, and extend it to equilibrium n-type and p-type doping.

This is an open computational investigation for a four-student course, not a fixed programming exercise. You have access to the full `physkit` library and may use any relevant models, solvers, unit systems, numerical tools, or visualization components. The intellectual work lies in choosing and connecting appropriate capabilities, checking their assumptions, validating the results, extending them where necessary, and explaining the resulting physical model.

The project ends at **equilibrium carrier statistics**. Transport is outside the scope: do not calculate mobility, conductivity, resistivity, Hall coefficients, Seebeck coefficients, diffusion coefficients, current density, or device performance.

The central modeling chain is

$$
H(\mathbf{k})
\longrightarrow E_n(\mathbf{k})
\longrightarrow \mathbf{M}^{*-1}
\longrightarrow g(E)
\longrightarrow f(E,\mu,T)
\longrightarrow n(T),p(T).
$$

Here $H(\mathbf{k})$ is a Bloch Hamiltonian, $E_n(\mathbf{k})$ is band $n$, $\mathbf{M}^*$ is an effective-mass tensor, $g(E)$ is the density of states, $f(E,\mu,T)$ is the Fermi–Dirac occupation, $\mu$ is the electronic chemical potential, and $n$ and $p$ are the equilibrium concentrations of electrons and holes.

At finite temperature, use **chemical potential** $\mu(T)$ or **Fermi level** rather than using “Fermi energy” indiscriminately. The Fermi energy is strictly a zero-temperature concept, although semiconductor literature often uses “Fermi level” for $\mu(T)$.

## Learning outcomes

After completing the assessment, you should be able to:

1. extract a band-edge effective-mass tensor from numerical band energies;
2. distinguish a full band model from its local parabolic effective-mass approximation;
3. distinguish a two-dimensional validation model from a three-dimensional bulk-semiconductor model;
4. construct conduction- and valence-band densities of states with correct dimensionality, units, spin degeneracy, and valley degeneracy;
5. calculate Fermi–Dirac occupations and solve charge-neutrality equations numerically;
6. calculate $\mu_i(T)$, $n_i(T)$, $n(T)$, and $p(T)$ without using transport equations;
7. reconstruct a defensible material model from tabulated data and published band diagrams;
8. identify when isotropic, parabolic, nondegenerate approximations are valid and when they fail;
9. implement incomplete donor and acceptor ionization using a literature-supported model; and
10. separate physical assumptions, numerical approximations, parameter uncertainty, and implementation error.

## Assessment conditions

- This is an **individual, open-note, open-computing-environment assessment**.
- You may use the full public [`physkit`](https://github.com/eragasa/projectkoios-physkit) codebase and the assessment starter repository. You are not expected to reimplement maintained package capabilities merely to demonstrate that you can reproduce them.
- Identify which capabilities are supplied by `physkit`, which configurations and analyses are your own, and what you changed or extended. Your scientific choices, material reconstruction, validation, interpretation, and defense must be your own.
- Students may discuss methods and share references. Any shared code or data must be identified; each submitted computational argument and interpretation must remain individually defensible.
- Your submitted computational record must run from a clean environment from beginning to end.
- All numerical quantities must carry documented units. Internally consistent SI units are recommended; final carrier concentrations must be reported in $\mathrm{cm^{-3}}$.
- Every material parameter must have a source, temperature, crystal phase or polytype, and definition.
- A parameter copied from a webpage without checking its cited reference is not automatically authoritative.
- The selected semiconductor must be a **three-dimensional bulk material**. Thin films, quantum wells, nanowires, and surface states are outside the common core.
- Because this is a small-cohort course, each student may pursue a different material, numerical strategy, or level of model complexity. The project declaration will be used to agree on an individually appropriate scope while preserving the common physical objectives.
- The sections below describe the evidence that should emerge from the investigation. They are not intended as a point-by-point recipe; equivalent methods and well-justified alternatives are welcome.

## Required project declaration

Before beginning the calculation, submit a one-page declaration containing:

| Item | Required statement |
|---|---|
| Material identity | Chemical formula and crystal phase or polytype |
| Sample model | Ideal bulk crystal, with the modeled temperature range |
| Gap model | Direct or indirect; $E_g(T)$ equation or tabulated values |
| Conduction edge | Location, number of equivalent valleys, mass model |
| Valence edge | Relevant bands, degeneracies, and mass model |
| Statistics | Full Fermi–Dirac integration, Boltzmann approximation, or both |
| n-type model | Dopant or donor-like defect and charge-state model |
| p-type model | Dopant or acceptor-like defect and charge-state model |
| Exclusions | Physics known to be outside the submitted model |

“SiC” is not a complete material identity. You must specify, for example, 3C-SiC, 4H-SiC, or 6H-SiC. “Bismuth telluride” is also insufficient unless you specify bulk crystalline $\mathrm{Bi_2Te_3}$ and state whether surface states are excluded.

---

## Part I: Validate the effective-mass reduction using `Bloch2D`

### I.1 Establish the parent band calculation

Use appropriate `physkit` `Bloch2D` capabilities to establish a trustworthy parent calculation. Explain enough of the direct and reciprocal lattice, boundary conditions, numerical representation, Hamiltonian, and eigensolution that another student could understand and reproduce the model. Include a meaningful Hermiticity or operator-consistency check.

Use both a high-symmetry path for interpretation and a two-dimensional Brillouin-zone sampling for extrema and state counting. A band path has zero area and cannot by itself support a DOS calculation. Clearly distinguish package-supplied behavior from your own configuration, analysis, and extensions.

### I.2 Locate and classify the band edges

Numerically locate one conduction-band minimum and one valence-band maximum. Do not assume that an extremum lies exactly at a plotted high-symmetry point. Report:

- the band index;
- the numerical $\mathbf{k}$ location;
- the edge energy;
- whether the extremum is isolated or symmetry-equivalent to other extrema; and
- the energy window used for the local fit.

Use a local quadratic model. For a conduction minimum at $\mathbf{k}_c$,

$$
E_c(\mathbf{k}_c+\mathbf{q})
\approx E_C+\frac{\hbar^2}{2}\mathbf{q}^{\mathsf T}
\mathbf{M}_c^{-1}\mathbf{q}.
$$

For a valence maximum at $\mathbf{k}_v$,

$$
E_v(\mathbf{k}_v+\mathbf{q})
\approx E_V-\frac{\hbar^2}{2}\mathbf{q}^{\mathsf T}
\mathbf{M}_v^{-1}\mathbf{q},
$$

where $\mathbf{M}_v$ is the positive hole-mass tensor. Extract the inverse-mass tensor from the Hessian,

$$
\left(\mathbf{M}^{-1}\right)_{ij}
=\frac{1}{\hbar^2}
\frac{\partial^2E}{\partial k_i\partial k_j}.
$$

Diagonalize the fitted tensor and report its principal masses and principal directions. A negative eigenvalue at a claimed conduction minimum, or a positive electronic curvature at a claimed valence maximum before conversion to hole mass, must be investigated.

### I.3 Validate the curvature extraction

Repeat the fit using at least three fitting radii or energy windows. Show the convergence or drift of the principal effective masses. Separate:

- error caused by the $\mathbf{k}$-mesh spacing;
- error caused by the quadratic-fit window; and
- physical nonparabolicity of the band.

For at least two directions through the extremum, overlay the numerical band energies and the fitted parabola. Report a quantitative residual, such as the maximum absolute energy error or root-mean-square error within the fit window.

### I.4 Validate density-of-states state counting

Calculate the two-dimensional density of states directly from the uniform $\mathbf{k}$ mesh using a normalized histogram, Gaussian broadening, or another documented kernel. Then calculate the analytic two-dimensional parabolic-band density of states implied by the fitted masses.

For an anisotropic two-dimensional valley with principal masses $m_1$ and $m_2$, the per-valley density-of-states mass is

$$
m_{d,2D}=\sqrt{m_1m_2}.
$$

With spin degeneracy $g_s$ and valley degeneracy $g_v$, the ideal parabolic conduction-band DOS per unit area is

$$
g_{c,2D}(E)=
\frac{g_sg_vm_{d,2D}}{2\pi\hbar^2}\,\Theta(E-E_C).
$$

Compare the numerical and analytic DOS near the selected edge. Explain the effects of finite mesh spacing, broadening, finite bandwidth, additional valleys, and nonparabolicity.

> [!CALLOUT]
>
> **Model boundary:** `Bloch2D` validates the numerical sequence from band curvature to DOS and occupation. It does not supply the three-dimensional DOS of the real bulk semiconductor used in Parts II–IV.

### Required outputs for Part I

- band structure on a high-symmetry path;
- band-edge map on a uniform $\mathbf{k}$ mesh;
- table of fitted tensor components and principal masses;
- fit-window convergence plot;
- numerical versus analytic 2D DOS plot; and
- a short conclusion stating the energy range over which the parabolic reduction is valid.

---

## Part II: Construct and validate an intrinsic three-dimensional effective-mass model

Use any suitable `physkit` capabilities for units, densities of states, Fermi–Dirac statistics, quadrature, root solving, data structures, or visualization. You may compose existing components, extend them, or implement a focused alternative where the library does not match your model. In every case, make the physical assumptions and validation evidence visible.

### II.1 Define the physical state model

Implement a bulk semiconductor with a conduction edge $E_C$, valence edge $E_V$, and band gap

$$
E_g=E_C-E_V.
$$

Use $E_V=0$ as the default internal energy reference unless another reference is explicitly justified. Your code must not depend on this arbitrary zero: shifting all energies and $\mu$ by a constant must leave $n$ and $p$ unchanged.

Represent each relevant band-edge valley by

$$
E_c(\mathbf{q})=E_C+
\frac{\hbar^2}{2}\mathbf{q}^{\mathsf T}\mathbf{M}_c^{-1}\mathbf{q},
$$

or

$$
E_v(\mathbf{q})=E_V-
\frac{\hbar^2}{2}\mathbf{q}^{\mathsf T}\mathbf{M}_v^{-1}\mathbf{q}.
$$

For a three-dimensional ellipsoid with principal masses $m_1,m_2,m_3$, define the per-valley density-of-states mass by

$$
m_d=(m_1m_2m_3)^{1/3}.
$$

If $g_v$ equivalent valleys are absorbed into a single quoted DOS mass, then

$$
m_{d,\mathrm{total}}=g_v^{2/3}m_d.
$$

Do not multiply by $g_v$ again after using a mass that already includes valley degeneracy.

### II.2 Implement the three-dimensional DOS

For one set of equivalent parabolic conduction valleys,

$$
g_c(E)=
\frac{g_sg_v}{4\pi^2}
\left(\frac{2}{\hbar^2}\right)^{3/2}
\sqrt{\det\mathbf{M}_c}\,
\sqrt{E-E_C}\,\Theta(E-E_C).
$$

For the valence band,

$$
g_v(E)=
\frac{g_sg_{v,h}}{4\pi^2}
\left(\frac{2}{\hbar^2}\right)^{3/2}
\sqrt{\det\mathbf{M}_v}\,
\sqrt{E_V-E}\,\Theta(E_V-E).
$$

If multiple inequivalent valleys or valence bands contribute, sum their DOS contributions explicitly. For degenerate heavy- and light-hole bands with isotropic masses,

$$
m_{v,\mathrm{DOS}}=
\left(m_{hh}^{3/2}+m_{lh}^{3/2}\right)^{2/3},
$$

provided the two maxima have the same energy and the spin factors are treated consistently.

Verify the units of $g(E)$ and show that integrating it produces a concentration per unit volume.

### II.3 Implement occupation and carrier concentration

Use the Fermi–Dirac function

$$
f(E,\mu,T)=\frac{1}{1+\exp[(E-\mu)/(k_BT)]}.
$$

Implement numerically stable evaluation for large positive or negative arguments. Then calculate

$$
n(\mu,T)=\int_{E_C}^{\infty}g_c(E)f(E,\mu,T)\,dE,
$$

and

$$
p(\mu,T)=\int_{-\infty}^{E_V}g_v(E)[1-f(E,\mu,T)]\,dE.
$$

The integration limits must be chosen by a convergence test rather than by an undocumented constant. Your code must support multiple valleys or bands by summing their contributions.

### II.4 Solve the intrinsic neutrality problem

For an ideal intrinsic semiconductor, determine $\mu_i(T)$ from

$$
n(\mu_i,T)-p(\mu_i,T)=0.
$$

Use a bracketing root solver. Document the bracket, tolerance, and residual. The result is unacceptable if the returned root does not satisfy charge neutrality to a stated tolerance.

Implement the nondegenerate analytic limit as an independent validation oracle:

$$
N_C(T)=2g_{v,e}
\left(\frac{2\pi m_{c,d}k_BT}{h^2}\right)^{3/2},
$$

$$
N_V(T)=2g_{v,h}
\left(\frac{2\pi m_{v,d}k_BT}{h^2}\right)^{3/2},
$$

$$
n_i(T)=\sqrt{N_CN_V}\exp\left[-\frac{E_g}{2k_BT}\right],
$$

and

$$
\mu_i(T)=\frac{E_C+E_V}{2}
+\frac{k_BT}{2}\ln\left(\frac{N_V}{N_C}\right).
$$

Compare the numerical Fermi–Dirac integrals with the nondegenerate formulas. Explain when the two agree and when the approximation fails.

### II.5 Required silicon control case

Before analyzing the selected material, run a silicon control calculation at $T=300\,\mathrm{K}$. Use a single internally consistent parameter source. Your calculation should recover approximately

$$
E_g\approx1.12\,\mathrm{eV},
\qquad
n_i\approx1.5\times10^{10}\,\mathrm{cm^{-3}},
$$

with $\mu_i$ close to, but not necessarily exactly at, midgap. Published values differ because of temperature, gap model, and effective-DOS conventions. A discrepancy must be traced to the input definitions rather than hidden by changing constants.

### II.6 Validation evidence

Build a compact set of automated or reproducible checks appropriate to your chosen implementation. Taken together, the evidence should address monotonic behavior of $n$ and $p$ with chemical potential, charge-neutrality residuals, invariance under a rigid energy shift, isotropic and anisotropic mass limits, the expected $T^{3/2}$ scaling of effective densities for constant masses, and convergence of numerical integration. You may use `physkit` tests and validation objects as part of this evidence, but explain what each check establishes for your model.

### Required outputs for Part II

- plot of $g_c(E)$ and $g_v(E)$ with $f(E,\mu_i,T)$ at one temperature;
- plot of $\mu_i(T)$ relative to $E_V$, midgap, and $E_C$;
- log-scale plot of $n_i(T)$;
- table comparing numerical and nondegenerate results at selected temperatures;
- silicon control-case table; and
- test results or a validation table.

---

## Part III: Reconstruct and analyze a real semiconductor from literature

### III.1 Source requirements

Choose one semiconductor and one specific bulk crystal structure. Assemble a sufficiently varied source set to reconstruct and challenge the model rather than copying a single table. A typical source set will include a curated database or handbook, primary and review literature, and a published band diagram or numerical band-data source. The exact mix may differ by material and should be justified in the project declaration.

The [Ioffe Institute semiconductor archive](https://www.ioffe.ru/SVA/NSM/Semicond/) is an acceptable starting point, not a substitute for checking its cited literature. If the material is not covered there, use an equivalent handbook, review, primary paper, or openly documented computational database.

Prepare a provenance table containing:

| Quantity | Symbol | Value or equation | Units | Temperature | Phase/polytype | Definition | Source |
|---|---|---:|---|---:|---|---|---|
| Band gap | $E_g$ |  | eV |  |  | direct/indirect |  |
| Electron mass data | $\mathbf{M}_c$ or $m_{c,d}$ |  | $m_0$ |  |  | curvature/DOS/other |  |
| Hole mass data | $\mathbf{M}_v$ or $m_{v,d}$ |  | $m_0$ |  |  | HH/LH/DOS/other |  |
| Valley degeneracy | $g_v$ |  | 1 |  |  | included in mass? |  |
| Conduction DOS | $N_C$ |  | $\mathrm{cm^{-3}}$ |  |  | measured/derived |  |
| Valence DOS | $N_V$ |  | $\mathrm{cm^{-3}}$ |  |  | measured/derived |  |
| Published intrinsic concentration | $n_i$ |  | $\mathrm{cm^{-3}}$ |  |  | method |  |

Do not mix masses, gaps, or valley counts from different crystal phases without constructing an explicit comparison.

### III.2 Use the published band diagram correctly

From a published band diagram, identify:

- $E_V$ and $E_C$;
- direct or indirect character of the gap;
- the location and multiplicity of relevant extrema;
- nearby bands that may contribute within the chosen temperature range; and
- evidence for anisotropy or nonparabolicity.

If numerical band data are unavailable, digitize local band-edge points from the figure and fit their curvature. Record the reciprocal-space axis, the conversion from reduced coordinates to $\mathrm{m^{-1}}$, the selected points, and the fit interval.

A single curvature extracted along one high-symmetry line does **not** determine a three-dimensional DOS mass. To calculate $n$ and $p$, you must obtain the remaining principal curvatures from literature, use a published DOS mass, or explicitly justify a symmetry relation. A high-symmetry path may also miss off-symmetry extrema. These limitations must be discussed.

### III.3 Construct the material model

Implement $E_g(T)$ when supported by literature. Recalculate, over a justified temperature range:

- $N_C(T)$ and $N_V(T)$;
- $g_c(E,T)$ and $g_v(E,T)$;
- $\mu_i(T)$;
- $n(T)$ and $p(T)$; and
- $n_i(T)$.

For the intrinsic model, verify numerically that $n=p$ at every reported temperature. Do not calculate $n$ from one approximation and $p$ from a different approximation unless the mismatch is itself the subject of the analysis.

### III.4 Compare against published values

Compare your model with at least one published value or curve for $N_C$, $N_V$, $\mu_i$, or $n_i$. Report the signed and percentage differences when meaningful. Discuss which discrepancies arise from:

- the band-gap model;
- mass definition;
- valley degeneracy;
- temperature dependence of mass;
- nonparabolicity;
- additional valence or conduction bands;
- digitization uncertainty;
- Fermi–Dirac versus Boltzmann statistics; and
- sample nonstoichiometry or unintentional doping in the experiment.

Perform a sensitivity analysis. At minimum, vary $E_g$, the electron DOS mass, and the hole DOS mass within literature ranges or reasonable uncertainty bounds and show their effects on $n_i$ and $\mu_i$.

## Material-specific pathways

### Path A: silicon carbide

Silicon carbide is allowed only when a single polytype is specified. **4H-SiC is the recommended SiC case** because useful band, DOS, and dopant data are available, but 3C-SiC or 6H-SiC may be used with consistent sources.

A defensible SiC investigation will normally address:

- an indirect gap;
- the correct conduction-valley multiplicity;
- an anisotropic electron mass or a properly defined DOS mass;
- a justified valence-band DOS mass;
- $E_g(T)$; and
- log-domain evaluation because $n_i$ can be extremely small near room temperature.

As a 4H-SiC sanity check, the Ioffe compilation gives approximately $E_g=3.23\,\mathrm{eV}$, $N_C=1.7\times10^{19}\,\mathrm{cm^{-3}}$, and $N_V=2.5\times10^{19}\,\mathrm{cm^{-3}}$ at $300\,\mathrm{K}$. These imply an intrinsic concentration of order $10^{-8}\,\mathrm{cm^{-3}}$, not a silicon-like value. Your exact result must follow your chosen source set.

Suitable doping studies include nitrogen or phosphorus donors and aluminum or boron acceptors. Incomplete ionization is essential for many SiC dopants; a full-ionization model at all temperatures is not acceptable. Site-dependent donor energies in hexagonal polytypes must be acknowledged when relevant.

### Path B: cadmium telluride

Bulk zinc-blende CdTe is the cleanest of the three example choices for the common effective-mass framework. It has a direct gap at $\Gamma$ and a relatively simple conduction-band edge, but its valence edge requires at least heavy- and light-hole consideration or a defensible combined DOS mass.

A defensible CdTe investigation will normally address:

- the direct $\Gamma$ gap and $E_g(T)$;
- one conduction-band edge with a literature-supported electron mass;
- heavy- and light-hole contributions or a published valence DOS mass;
- comparison with the published intrinsic $n_i(T)$ or $\mu_i(T)$ literature; and
- a discussion of native-defect compensation.

The study by Ching-Hua Su, *Journal of Applied Physics* **103**, 084903 (2008), DOI [10.1063/1.2899087](https://doi.org/10.1063/1.2899087), provides a particularly suitable comparison because it analyzes the band gap, intrinsic carrier concentration, and Fermi level of bulk CdTe from $304$ to $1067\,\mathrm{K}$ using a Kane conduction band and light- and heavy-hole valence bands.

Suitable doping studies include a chlorine or indium donor model and an arsenic, phosphorus, or copper acceptor model. CdTe doping is not adequately described by declaring a nominal dopant concentration and assuming one free carrier per dopant: compensation by vacancies, antisites, amphoteric behavior, and defect complexes may be important. Your selected model may remain a one-level equilibrium model, but its limitations must be compared with the defect literature.

### Path C: bulk bismuth telluride

Bulk crystalline $\mathrm{Bi_2Te_3}$ is an **advanced model-validity case**, not an ordinary isotropic parabolic-band exercise. It is a narrow-gap, layered, strongly anisotropic material with strong spin–orbit coupling, off-symmetry extrema, multivalley band edges, and topological surface states. In this assessment, model only three-dimensional **bulk** carriers and explicitly exclude surface-state transport and two-dimensional Dirac carriers.

A defensible $\mathrm{Bi_2Te_3}$ investigation will normally address:

- an SOC-consistent bulk gap;
- the relevant bulk-valley multiplicity;
- anisotropic DOS masses or a published full bulk DOS;
- a test of parabolicity over the thermal energy window;
- full Fermi–Dirac statistics when $\mu$ approaches or enters a band; and
- a clear separation between the ideal intrinsic bulk model and commonly measured, unintentionally doped samples.

A six-valley model is often required for bulk carrier statistics. The band extrema need not lie on the simplest high-symmetry path, so a plotted path alone cannot establish the DOS. Strong spin–orbit coupling substantially changes the gap and curvature; parameters from a calculation without SOC cannot be mixed with an SOC band gap. A useful modern overview is Witting *et al.*, “The Thermoelectric Properties of Bismuth Telluride,” *Advanced Electronic Materials* **5**, 1800904 (2019), DOI [10.1002/aelm.201800904](https://doi.org/10.1002/aelm.201800904). The band-edge uncertainty and off-symmetry extrema should be treated as part of the scientific analysis, not hidden.

For doping, choose one literature-supported donor-like and one acceptor-like substitution, vacancy, or antisite model. Because real $\mathrm{Bi_2Te_3}$ samples are frequently nonstoichiometric and degenerate, a defect-chemistry model may be more defensible than a hydrogenic shallow-dopant model. You are not required to calculate defect formation energies, but you must state what fixes the modeled defect concentration and which charge states enter neutrality.

---

## Part IV: Extend the material model to equilibrium n-type and p-type doping

### IV.1 Select literature-supported dopants

Choose one donor-like species and one acceptor-like species for the selected material. For each, report:

- chemical identity and lattice site or defect type;
- energy level relative to $E_C$ or $E_V$;
- charge states included;
- degeneracy factor and its convention;
- concentration range;
- experimental or theoretical source; and
- known compensation or amphoteric behavior.

Do not use the hydrogenic donor formula as a source of the dopant ionization energy when a measured material-specific value is available. You may use it as a comparison and discuss why it succeeds or fails.

### IV.2 Implement incomplete ionization

For a single donor level $E_D$ and acceptor level $E_A$, one common convention is

$$
N_D^+(\mu,T)=
\frac{N_D}{1+g_D\exp[(\mu-E_D)/(k_BT)]},
$$

$$
N_A^-(\mu,T)=
\frac{N_A}{1+g_A\exp[(E_A-\mu)/(k_BT)]}.
$$

Degeneracy conventions differ across references. Define $g_D$ and $g_A$ exactly as used in your equations and test the low- and high-temperature limits. Copying a degeneracy factor without its convention is an error.

For simultaneous donors and acceptors, solve

$$
p(\mu,T)+N_D^+(\mu,T)
-n(\mu,T)-N_A^-(\mu,T)=0.
$$

Your solver must return $\mu$, $n$, $p$, $N_D^+$, and $N_A^-$ and must report the neutrality residual.

### IV.3 Representative doping cases

Choose a set of cases broad enough to demonstrate intrinsic, n-type, p-type, and compensation behavior. A natural set includes:

1. intrinsic material: $N_D=N_A=0$;
2. uncompensated n-type material;
3. compensated n-type material;
4. uncompensated p-type material; and
5. compensated p-type material.

If the material's documented defect chemistry makes another comparison more informative, propose it in the project declaration.

Select concentration and temperature ranges from the material literature. The range must expose, where physically present:

- freeze-out;
- an extrinsic regime;
- onset of intrinsic excitation; and
- degeneracy or breakdown of the Boltzmann approximation.

If a regime is not reached within the material’s stable or documented temperature range, say so rather than extending the model without justification.

### IV.4 Required analysis

Produce:

- $\mu(T)$ relative to $E_V$, $E_C$, $E_D$, and $E_A$;
- $n(T)$ and $p(T)$ on a logarithmic scale;
- donor and acceptor ionization fractions versus temperature;
- $n$ and $p$ versus nominal dopant concentration at one or more temperatures;
- neutrality residuals over the complete sweep; and
- equilibrium band diagrams for representative intrinsic, n-type, and p-type cases.

The band diagrams must show energy, occupation context, and charge state; they must not imply spatial band bending because this project contains no Poisson equation or device geometry.

Compare the incomplete-ionization model with the fully ionized approximation. Identify the parameter range over which the simpler model gives less than 5% error in the majority-carrier concentration. For a narrow-gap or degenerately doped case, also compare Fermi–Dirac and Boltzmann results.

### IV.5 Scientific conclusion

End with a defensible answer to the following question:

> Over what temperature and concentration range is the chosen effective-mass plus equilibrium dopant model an adequate description of the selected semiconductor, and which omitted physics is most likely to invalidate it first?

Your answer must cite numerical evidence from your calculations and evidence from the literature.

---

## Computational organization and use of `physkit`

Organize the work so that the physical model, numerical choices, calculated results, validation evidence, and plots can be distinguished. You may use `physkit` objects directly, compose them with small student-authored functions, or develop a more specialized model where appropriate. There is no required class hierarchy or function naming scheme.

Avoid hiding independent physics calculations inside plotting cells. Preserve enough parameter metadata and provenance that a result can be traced back to its material assumptions, package capabilities, and numerical settings. Prefer the simplest organization that remains readable, reproducible, and extensible for your chosen investigation.

## Suggested submission structure

```text
semiconductor-assessment/
├── README.md
├── environment.yml or requirements.txt
├── data/
│   ├── material-parameters.yaml
│   ├── digitized-band-data.csv
│   └── source-notes.md
├── notebooks/
│   ├── 01_bloch2d_validation.ipynb
│   ├── 02_intrinsic_effmass.ipynb
│   ├── 03_material_analysis.ipynb
│   └── 04_doping_analysis.ipynb
├── src/
│   └── semiconductor_model.py
├── tests/
│   └── test_semiconductor_model.py
└── report.pdf
```

## Submission package

Submit:

1. a report of approximately 3,000–5,000 words, excluding equations, captions, tables, and references;
2. one or more executable notebooks that present a clearly organized computational argument;
3. any student-authored reusable Python source code used by the analysis;
4. automated tests, validation tables, or equivalent reproducible checks;
5. the material-parameter provenance table;
6. digitized band data, if used;
7. all plots in vector or high-resolution form; and
8. a short oral defense.

The report must distinguish:

- **physics:** the band and carrier-statistics assumptions;
- **mathematics:** DOS integrals and neutrality equations;
- **numerics:** meshes, fits, quadrature, root solving, and convergence;
- **implementation:** data structures, functions, tests, and units; and
- **evidence:** comparison with independent published values.

## Oral defense

Prepare a 10-minute defense followed by questions. Your defense must contain one compact model diagram, one validation table, one decisive plot, and one conclusion. Be prepared to answer:

1. Why can a high-symmetry band path not be used directly as a DOS integration mesh?
2. How is a curvature mass different from a DOS effective mass?
3. Where did spin and valley degeneracy enter your calculation?
4. Why is intrinsic $\mu_i$ not always exactly at midgap?
5. What numerical evidence shows that your neutrality solver succeeded?
6. When does the Boltzmann approximation fail?
7. Which parameter most strongly controls $n_i$?
8. What part of the `Bloch2D` validation transfers to the real material, and what part does not?
9. Why does the chosen dopant not necessarily produce one free carrier per dopant atom?
10. Which omitted physical effect would invalidate your model first?

## Evaluation approach

This assessment is evaluated holistically rather than through a point-by-point rubric. The four parts form one connected scientific argument, and different materials or modeling choices may require different kinds of evidence.

The evaluation will consider the overall quality of:

- the physical model and the consistency of its assumptions;
- the connection between equations, numerical methods, and implemented code;
- reproducibility, units, provenance, and transparent state counting;
- validation against limiting cases, convergence studies, and independent data;
- interpretation of discrepancies, uncertainty, and model limitations;
- the clarity of the written report, notebook organization, and oral defense; and
- the student's ability to explain and defend the choices made in the project.

A successful submission should present a coherent path from the `Bloch2D` calculation to a three-dimensional intrinsic model and then to equilibrium donor and acceptor cases. It should provide enough numerical and documentary evidence to establish that the calculations are reproducible, charge neutral where required, and appropriate for the selected material. Alternative formulations are acceptable when their assumptions are explicit and their validity is demonstrated.

Depth, scope, and suitable validation evidence will be discussed through the project declaration and short progress consultations. More complex materials are not expected to produce the same model form as simpler materials; they are expected to show an appropriate level of judgment about what can and cannot be represented reliably.

## Starting references

1. Ioffe Institute, [New Semiconductor Materials: Characteristics and Properties](https://www.ioffe.ru/SVA/NSM/Semicond/).
2. Ioffe Institute, [Silicon carbide band structure and carrier concentration](https://www.ioffe.ru/SVA/NSM/Semicond/SiC/bandstr.html).
3. Ching-Hua Su, “Energy band gap, intrinsic carrier concentration, and Fermi level of CdTe bulk crystal between 304 and 1067 K,” *Journal of Applied Physics* **103**, 084903 (2008), [DOI: 10.1063/1.2899087](https://doi.org/10.1063/1.2899087).
4. D. T. F. Marple, “Effective Electron Mass in CdTe,” *Physical Review* **129**, 2466 (1963), [DOI: 10.1103/PhysRev.129.2466](https://doi.org/10.1103/PhysRev.129.2466).
5. I. G. Ivanov, A. Henry, and E. Janzén, “Ionization energies of phosphorus and nitrogen donors and aluminum acceptors in 4H-SiC from donor-acceptor pair emission,” *Physical Review B* **71**, 241201(R) (2005), [DOI: 10.1103/PhysRevB.71.241201](https://doi.org/10.1103/PhysRevB.71.241201).
6. I. T. Witting *et al.*, “The Thermoelectric Properties of Bismuth Telluride,” *Advanced Electronic Materials* **5**, 1800904 (2019), [DOI: 10.1002/aelm.201800904](https://doi.org/10.1002/aelm.201800904).
7. I. A. Nechaev and E. V. Chulkov, “Quasiparticle band gap in the topological insulator Bi$_2$Te$_3$,” *Physical Review B* **88**, 165135 (2013), [DOI: 10.1103/PhysRevB.88.165135](https://doi.org/10.1103/PhysRevB.88.165135).

These references are starting points. The submitted material model must use the original sources most appropriate to its selected phase, temperature range, mass definition, and dopant species.
