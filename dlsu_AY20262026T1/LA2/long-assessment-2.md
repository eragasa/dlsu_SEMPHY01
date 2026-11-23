# Long Assessment 2: From Equilibrium Carrier Statistics to Transport and p–n Junction Prediction

## Purpose and scope

In Long Assessment 1, you constructed a material-specific equilibrium model connecting band structure, effective masses, density of states, Fermi–Dirac occupation, charge neutrality, and incomplete dopant ionization. In Long Assessment 2, you will use that model as the initial state of a transport and device calculation.

The central question is:

> Given the band-edge, carrier-statistics, and dopant results established in LA1, what low-field transport and p–n-junction behavior should the selected semiconductor exhibit, and where does the predictive chain fail because additional physics or data are required?

The assessment connects

$$
\begin{aligned}
&\text{LA1 band and equilibrium model}\\
&\quad\longrightarrow
\{E_g,\mathbf M^*,n_i,\mu,n,p,N_D^+,N_A^-\}\\
&\quad\longrightarrow
\{\mu_n,\mu_p,D_n,D_p,\sigma,\rho,R_H\}\\
&\quad\longrightarrow
\{\tau_n,\tau_p,L_n,L_p,U\}\\
&\quad\longrightarrow
\{V_{\mathrm{bi}},W,E_{\max},C_j\}\\
&\quad\longrightarrow
\{J(V),C_d,Q_s,t_{\mathrm{rr}},\text{operating limits}\}.
\end{aligned}
$$

Here $\mathbf M^*$ is an effective-mass tensor; $n_i$ is intrinsic carrier concentration; $\mu$ is chemical potential; $n$ and $p$ are electron and hole densities; $N_D^+$ and $N_A^-$ are ionized donor and acceptor densities; $\mu_n$ and $\mu_p$ are drift mobilities; $D_n$ and $D_p$ are diffusivities; $\sigma$ is conductivity; $\rho$ is resistivity; $R_H$ is Hall coefficient; $\tau_n$ and $\tau_p$ are minority-carrier lifetimes; $L_n$ and $L_p$ are diffusion lengths; $U$ is net recombination rate; $V_{\mathrm{bi}}$ is built-in potential; $W$ is depletion width; $E_{\max}$ is peak electric field; $C_j$ is depletion capacitance; $J(V)$ is junction current density; $C_d$ is diffusion capacitance; $Q_s$ is stored minority charge; and $t_{\mathrm{rr}}$ is a reverse-recovery time.

This chain contains two distinct kinds of prediction:

1. **State propagation:** LA1 results such as $n$, $p$, $n_i$, dopant ionization, band gap, and effective mass are direct inputs to transport and junction equations.
2. **New constitutive physics:** mobility, scattering time, recombination lifetime, capture cross section, surface recombination velocity, impact-ionization coefficient, and contact resistance cannot be obtained from LA1 carrier statistics alone. They require new measurements, literature models, or explicitly declared assumptions.

A successful submission must preserve this distinction. Do not infer an absolute mobility from effective mass alone. Do not infer a lifetime from carrier concentration alone. Do not claim a predictive breakdown voltage without a defensible high-field ionization or tunneling model.

## Learning outcomes

After completing the assessment, you should be able to:

1. transfer a validated equilibrium material model into a transport calculation without silently changing its conventions;
2. explain which transport quantities follow from LA1 and which require new constitutive data;
3. distinguish drift mobility, conductivity mobility, and Hall mobility;
4. calculate conductivity, resistivity, Hall coefficient, and diffusivity for one- and two-carrier regimes;
5. test the nondegenerate Einstein relation and identify when a generalized relation is required;
6. construct a documented recombination and minority-carrier lifetime model;
7. solve steady and transient one-dimensional carrier-continuity equations numerically;
8. calculate p–n-junction built-in potential from the LA1 equilibrium states;
9. derive depletion width, field, potential, and capacitance for an abrupt junction;
10. compare the depletion approximation with a nonlinear Poisson or Poisson–Boltzmann solution;
11. calculate ideal and nonideal diode current from material-specific equilibrium and transport inputs;
12. calculate stored charge, diffusion capacitance, and a switching or recovery response;
13. quantify uncertainty propagation from LA1 and newly introduced transport parameters; and
14. state the bias, temperature, concentration, geometry, and field range over which the resulting device prediction is defensible.

## Assessment conditions

- This is an **individual, open-note, open-computing-environment assessment**.
- Continue with the same exact bulk semiconductor phase or polytype studied in LA1 unless a change is approved in the project declaration.
- Use the submitted LA1 material model as an input artifact. Do not retype selected LA1 results into plotting cells.
- The LA2 computational record must identify the exact LA1 revision, file hashes, or archived submission from which the inputs were obtained.
- A correction to LA1 is allowed, but the correction, its reason, and its effect on all downstream LA2 predictions must be documented.
- Use standard NumPy, SciPy, pandas, Matplotlib, Jupyter, and pytest dependencies. The submission must not depend on Project Koios or PhysKit.
- Use SI units internally unless a different system is explicitly justified. Report semiconductor concentrations in $\mathrm{cm^{-3}}$ as well as $\mathrm{m^{-3}}$ where comparison with literature requires it.
- Every new material or device parameter must have a definition, unit, temperature, concentration range, crystallographic direction when relevant, and source.
- Every p–n-junction calculation must state the junction orientation, geometry, area, neutral-region lengths, donor density, acceptor density, and temperature.
- Numerical solutions must include convergence, residual, conservation, positivity, or limiting-case evidence appropriate to the equation being solved.
- The project predicts a **model device**, not an unspecified commercial diode. Device geometry and boundary conditions are part of the model and must be declared.

## Provided silicon and gallium-arsenide examples

Four executed notebooks under `notebooks/` demonstrate the computational
sequence using silicon and gallium arsenide:

1. `01_si_gaas_transport.ipynb` — mobility, conductivity, resistivity, Hall
   response, and diffusivity;
2. `02_si_gaas_recombination.ipynb` — diffusion length, steady minority
   profiles, and transient charge removal;
3. `03_si_gaas_pn_electrostatics.ipynb` — built-in potential, depletion
   electrostatics, capacitance, and C–V extraction; and
4. `04_si_gaas_diode_prediction.ipynb` — finite-base saturation current,
   ideal J–V response, stored charge, and diffusion capacitance.

The notebooks begin with direct NumPy/SciPy calculations and then compare them
with the readable local implementation in `notebooks/la2_models.py`. The Si and
GaAs parameter sets are instructional controls. They demonstrate method,
validation, units, and model tracing; they are not automatically valid sources
for a student's LA1 material or for every Si or GaAs sample. Students who used
Si or GaAs in LA1 must still replace the example parameters with values
consistent with their own phase, orientation, temperature, concentration,
sample condition, and cited literature.

### Lecture-note reading map

The common-core derivations and numerical methods are developed in Chapters
7--21 of the course notes. Chapter 12 treats mobility, conductivity, Hall
response, and classical and generalized Einstein relations. Chapter 13 treats
recombination and finite-difference continuity. Chapters 16--20 treat junction
electrostatics, C--V response, current, charge storage, switching, and
breakdown. Chapter 21 explicitly connects the LA1 equilibrium-state record to
built-in voltage, ionized depletion charge, degenerate diffusion, calibration,
validation, sensitivity, and uncertainty.

Material-specific mobility, lifetime, surface, contact, and high-field
parameters must still be obtained from independent literature.

## Required LA1-to-LA2 transfer package

Before beginning new transport calculations, create an immutable transfer package under `data/la1-transfer/` containing:

1. the LA1 material-parameter file;
2. the LA1 source/provenance notes;
3. the student-authored equilibrium-semiconductor module used for the final LA1 results;
4. a machine-readable state table covering the temperatures and doping cases needed in LA2;
5. the LA1 validation summary; and
6. a manifest containing file hashes or the LA1 repository commit.

The state table must contain at least:

| Column | Meaning | Required unit |
|---|---|---|
| `material` | Chemical formula and exact phase/polytype | text |
| `temperature_K` | Absolute temperature | K |
| `case_id` | Intrinsic, n-type, p-type, or compensated state identifier | text |
| `donor_m3` | Nominal donor density | $\mathrm{m^{-3}}$ |
| `acceptor_m3` | Nominal acceptor density | $\mathrm{m^{-3}}$ |
| `ionized_donor_m3` | $N_D^+$ | $\mathrm{m^{-3}}$ |
| `ionized_acceptor_m3` | $N_A^-$ | $\mathrm{m^{-3}}$ |
| `electron_m3` | $n$ | $\mathrm{m^{-3}}$ |
| `hole_m3` | $p$ | $\mathrm{m^{-3}}$ |
| `chemical_potential_eV` | $\mu$ in the LA1 energy reference | eV |
| `intrinsic_m3` | $n_i$ | $\mathrm{m^{-3}}$ |
| `band_gap_eV` | $E_g$ | eV |
| `neutrality_residual_m3` | Signed charge-neutrality residual divided by $q$ | $\mathrm{m^{-3}}$ |

Include mass tensors, valley degeneracies, band-edge energies, and dopant levels in the accompanying parameter file rather than duplicating them in every state-table row.

### Transfer audit

Your first notebook must verify that:

- all required columns and units are present;
- the transferred neutrality residual remains within the LA1 tolerance;
- $np=n_i^2$ in the nondegenerate equilibrium cases where the mass-action law applies;
- a rigid shift of all band energies and chemical potentials does not change $n$, $p$, or any LA2 prediction;
- the selected n-side and p-side states lie within the validated temperature and concentration range of LA1; and
- the source files match the recorded manifest.

If the transfer audit fails, repair the LA1 model before continuing.

---

## Required project declaration

Submit a two-page declaration before the main calculation.

| Item | Required statement |
|---|---|
| Material continuity | Exact material and phase/polytype inherited from LA1 |
| LA1 revision | Commit, archive identifier, or file hashes |
| Temperature range | Common range supported by equilibrium, transport, and recombination data |
| Junction design | Orientation, area, $N_A$, $N_D$, neutral-region lengths, and contact assumptions |
| Mobility model | Electron and hole model, directions, variables, and sources |
| Hall model | One-carrier or two-carrier treatment and Hall-factor convention |
| Recombination model | Direct, SRH, Auger, surface, or measured-lifetime model |
| Lifetime inputs | Values or equations, injection range, and sources |
| Electrostatics | Depletion approximation and numerical comparison to be used |
| Current model | Long-base, finite-base, depletion recombination, series resistance, or other corrections |
| Dynamic prediction | Stored charge, diffusion capacitance, transient diffusion, or reverse recovery |
| Operating-limit model | Breakdown, tunneling, velocity saturation, heating, or documented exclusion |
| Validation data | Independent transport or diode data selected for comparison |
| Exclusions | Physics not represented in the submitted model |

---

# Part I: Transfer the LA1 state and construct a homogeneous transport model

## I.1 Explain what transfers from LA1

Construct a model map showing the inputs inherited from LA1 and the new transport inputs introduced in LA2.

At minimum, explain how the following LA1 quantities enter LA2:

- $E_g(T)$ controls $n_i(T)$ and influences reverse leakage and high-field limits;
- effective masses and valley structure constrain transport anisotropy and density-of-states conventions;
- $n(T)$ and $p(T)$ enter conductivity and Hall response;
- $N_D^+(T)$ and $N_A^-(T)$ determine the fixed depletion charge under the adopted ionization approximation;
- the n-type and p-type chemical potentials determine the contact potential; and
- compensation changes majority-carrier density even at fixed nominal doping.

Also explain why the following do not follow from LA1 without new physics:

- scattering time and mobility;
- minority-carrier lifetime;
- capture cross sections and center density;
- surface recombination velocity;
- contact resistance;
- impact-ionization coefficients; and
- heat removal.

## I.2 Mobility and scattering model

Construct electron and hole mobility models over the common temperature and concentration range. You may use:

- a measured mobility table with documented interpolation;
- a literature empirical model;
- a Matthiessen-rule model combining lattice and ionized-impurity terms;
- a relaxation-time model with independently justified scattering times; or
- another documented constitutive relation.

For one principal direction,

$$
\mu_n=\frac{q\tau_{m,n}}{m_{c,n}^*},
\qquad
\mu_p=\frac{q\tau_{m,p}}{m_{c,p}^*},
$$

where $\tau_m$ is momentum-relaxation time. This equation does not determine $\mu$ until $\tau_m$ is supplied.

If the material is anisotropic, use mobility tensors or clearly identified principal-direction mobilities. Do not combine a DOS mass with a directional mobility formula without justification.

Compare your mobility model with at least one independent measured value or curve. Report signed and percentage differences at representative state points.

## I.3 Conductivity and resistivity

Calculate

$$
\sigma=q(n\mu_n+p\mu_p),
\qquad
\rho=\sigma^{-1}
$$

for the intrinsic, n-type, p-type, and compensated cases transferred from LA1. Identify the carrier and mechanism controlling each temperature regime.

Required analyses:

- conductivity and resistivity versus temperature;
- conductivity versus nominal dopant concentration;
- comparison of incomplete and complete ionization;
- sensitivity to mobility parameters; and
- separation of changes caused by carrier density from changes caused by mobility.

## I.4 Hall response

For a single dominant carrier,

$$
R_H=\frac{s r_H}{qN},
$$

where $s=-1$ for electrons and $s=+1$ for holes, $N$ is the dominant carrier density, and $r_H$ is the Hall scattering factor.

When both carriers contribute,

$$
R_H=rac{p\mu_p^2-n\mu_n^2}
{q(p\mu_p+n\mu_n)^2}
$$

for unit Hall factors. If nonunit electron and hole Hall factors are used, define the adopted convention and modify the expression consistently.

Calculate the apparent Hall density and Hall mobility. Identify any temperature or compensation range in which Hall sign, majority carrier, and net dopant type do not give the same simple interpretation.

## I.5 Diffusion and the Einstein relation

In the nondegenerate limit,

$$
D_n=\mu_nV_T,
\qquad
D_p=\mu_pV_T,
\qquad
V_T=\frac{k_BT}{q}.
$$

Test whether the selected states are nondegenerate before applying this relation. For a degenerate state, use a generalized Einstein relation or explicitly bound the error introduced by the classical approximation.

### Required outputs for Part I

- LA1 transfer manifest and audit table;
- inherited-versus-new-input model map;
- mobility and scattering provenance table;
- $\mu_n(T,N)$ and $\mu_p(T,N)$ plots;
- $\sigma(T)$ and $\rho(T)$ plots for representative LA1 states;
- Hall coefficient, apparent Hall density, and Hall mobility plots;
- Einstein-relation validity table; and
- uncertainty or sensitivity analysis separating carrier-density and mobility effects.

---

# Part II: Construct and validate the nonequilibrium carrier model

## II.1 Select a recombination model

Choose a model appropriate to the selected material and injection range. The total bulk rate may contain

$$
U_{\mathrm{bulk}}=U_{\mathrm{dir}}+U_{\mathrm{SRH}}+U_{\mathrm A},
$$

where direct band-to-band, Shockley–Read–Hall, and Auger terms are included only when supported by the material and operating conditions.

For a single SRH center,

$$
U_{\mathrm{SRH}}=
\frac{np-n_i^2}
{\tau_p(n+n_1)+\tau_n(p+p_1)}.
$$

Define every center energy, degeneracy, capture coefficient or lifetime, and density. State whether the parameters are measured effective values or derived microscopic quantities.

A measured minority lifetime may be used as the common core when a complete center model is unavailable. In that case, identify the injection level, sample condition, temperature, and surface treatment associated with the measurement.

## II.2 Low-level injection and diffusion length

For a p-type neutral region under low-level electron injection,

$$
\frac{\partial\Delta n}{\partial t}
=D_n\frac{\partial^2\Delta n}{\partial x^2}
-v_d\frac{\partial\Delta n}{\partial x}
-\frac{\Delta n}{\tau_n}+G(x,t),
$$

with

$$
L_n=\sqrt{D_n\tau_n}.
$$

Construct the corresponding hole model for the n-type region. Verify low-level injection against the LA1 majority-carrier densities rather than against nominal dopant density alone.

## II.3 Required numerical continuity calculation

Implement a one-dimensional finite-difference continuity solver. The common core must support:

- diffusion;
- first-order lifetime recombination or the selected nonlinear rate;
- a documented source or initial condition;
- Dirichlet, flux, or surface-recombination boundary conditions; and
- time evolution with an implicit method.

Drift is optional in the common core and required when used in the final physical interpretation.

Validate the solver against at least two of the following:

1. uniform exponential decay;
2. the steady $\exp(-x/L)$ minority profile;
3. the Haynes–Shockley drifting Gaussian packet;
4. finite-base hyperbolic-sine solution; or
5. a manufactured solution.

Required numerical checks:

- grid refinement;
- time-step refinement;
- nonnegative carrier density;
- conservation of generated, recombined, stored, and boundary-extracted charge; and
- residual evaluation for a steady solution.

## II.4 Surface recombination or finite-base boundary physics

Include at least one nonideal boundary. A surface-recombination boundary may be written

$$
-D\frac{\partial\Delta c}{\partial n}=S\Delta c,
$$

where $S$ is surface recombination velocity in $\mathrm{m\,s^{-1}}$. Alternatively, use a finite neutral region terminated by an equilibrium ohmic contact. Justify the boundary physically.

### Required outputs for Part II

- recombination and lifetime provenance table;
- injection-level validity plot;
- diffusion lengths versus temperature;
- analytic-versus-numerical continuity validation;
- transient carrier profiles;
- charge-conservation table;
- grid and time-step convergence evidence; and
- sensitivity of diffusion length or stored charge to lifetime uncertainty.

---

# Part III: Predict p–n-junction electrostatics and capacitance

## III.1 Define the model junction from LA1 states

Construct an abrupt homojunction using one n-type and one p-type equilibrium state from LA1 at the same temperature. State:

- nominal $N_D$ and $N_A$;
- calculated $N_D^+$ and $N_A^-$;
- equilibrium $n$ and $p$ on each side;
- chemical potentials in the common LA1 energy reference;
- area and orientation;
- neutral-region lengths; and
- whether complete ionization, nondegenerate statistics, and one-dimensional electrostatics are valid.

Do not use nominal dopant density as fixed depletion charge when LA1 predicts substantial incomplete ionization without analyzing the resulting error.

## III.2 Built-in potential from the transferred equilibrium states

For a homojunction whose isolated bulk states share the same band-edge reference, calculate

$$
qV_{\mathrm{bi}}=\mu_n^{\mathrm{bulk}}-\mu_p^{\mathrm{bulk}}.
$$

Independently calculate the nondegenerate approximation

$$
V_{\mathrm{bi}}
=V_T\ln\left(\frac{n_{n0}p_{p0}}{n_i^2}\right)
\approx V_T\ln\left(\frac{N_DN_A}{n_i^2}\right)
$$

only when the corresponding assumptions hold. Here $n_{n0}$ is the equilibrium electron density in the n-type bulk and $p_{p0}$ is the equilibrium hole density in the p-type bulk.

Compare:

1. the chemical-potential result;
2. the carrier-density result;
3. the fully ionized textbook result; and
4. a numerical Poisson or Poisson–Boltzmann result.

Explain every discrepancy.

## III.3 Abrupt-junction depletion solution

For fixed ionized dopant densities,

$$
N_A^-x_p=N_D^+x_n,
$$

and

$$
W=x_p+x_n
=\sqrt{\frac{2\varepsilon_s}{q}(V_{\mathrm{bi}}-V)
\left(\frac{1}{N_A^-}+\frac{1}{N_D^+}\right)}.
$$

Calculate $x_p$, $x_n$, $W$, charge density, electric field, peak field, and potential for equilibrium, one forward bias, and at least three reverse biases.

## III.4 Numerical electrostatics

Implement a finite-difference Poisson or Poisson–Boltzmann calculation. The numerical model must state:

- the charge-density expression;
- treatment of incomplete ionization;
- boundary conditions;
- grid and domain size;
- nonlinear solver and tolerances; and
- residual normalization.

Compare numerical and depletion-approximation potential, field, and integrated charge. Verify global charge neutrality and boundary-field behavior.

## III.5 Depletion capacitance and C–V prediction

Calculate

$$
C_j=\frac{\varepsilon_s}{W}
$$

per unit area and the total capacitance $A_jC_j$. For an approximately one-sided abrupt junction, use

$$
N(W)=\frac{2}{q\varepsilon_s}
\left[\frac{\mathrm d(1/C_j^2)}{\mathrm dV_R}\right]^{-1}
$$

as an inverse check. Determine whether C–V extraction recovers the ionized dopant density, nominal dopant density, or neither under your assumptions.

### Required outputs for Part III

- table of the two transferred LA1 bulk states;
- four-way built-in-potential comparison;
- charge, field, potential, and band diagrams at selected biases;
- numerical-versus-depletion comparison;
- nonlinear residual and charge-neutrality checks;
- $C_j(V_R)$ and $1/C_j^2$ plots;
- recovered doping profile; and
- uncertainty propagation from $n_i$, dopant ionization, and permittivity to $V_{\mathrm{bi}}$, $W$, and $C_j$.

---

# Part IV: Predict diode current, stored charge, and switching behavior

## IV.1 Minority-carrier boundary conditions

Under the nondegenerate law of the junction,

$$
\Delta n_p(0)=n_{p0}\left(e^{V/V_T}-1\right),
\qquad
\Delta p_n(0)=p_{n0}\left(e^{V/V_T}-1\right),
$$

where $n_{p0}$ and $p_{n0}$ must come from the LA1 equilibrium states. State the voltage and injection range over which these boundary conditions remain valid.

## IV.2 Ideal finite-base diode current

For neutral-region lengths $w_p$ and $w_n$ terminated by equilibrium contacts, calculate

$$
J_s=q\left[
\frac{D_nn_{p0}}{L_n}\coth\left(\frac{w_p}{L_n}\right)
+
\frac{D_pp_{n0}}{L_p}\coth\left(\frac{w_n}{L_p}\right)
\right]
$$

and

$$
J(V)=J_s\left(e^{V/V_T}-1\right).
$$

Compare the long-base and finite-base limits. Identify which side dominates $J_s$ and trace that dominance to LA1 carrier densities and LA2 transport parameters.

## IV.3 Nonideal current

Include at least two of the following, selected according to available evidence:

- depletion-region generation under reverse bias;
- depletion-region SRH recombination under forward bias;
- high-level injection;
- series resistance;
- surface leakage;
- band-to-band tunneling;
- temperature dependence; or
- a documented material-specific correction.

A common forward model is

$$
J=J_{01}\left(e^{V/V_T}-1\right)
+J_{02}\left(e^{V/(2V_T)}-1\right),
$$

possibly with junction voltage replaced by $V-JR_s$. Every fitted or sourced parameter must be identified.

## IV.4 Numerical minority-current validation

Use the Part II finite-difference solver to calculate the steady minority profiles on both sides of the junction. Extract the edge fluxes and compare them with the analytic finite-base current. Verify current continuity and convergence of the numerical gradient.

## IV.5 Stored charge and diffusion capacitance

Calculate the minority stored charge

$$
Q_s=q\int\Delta c(x)\,\mathrm dx
$$

and diffusion capacitance

$$
C_d=\frac{\mathrm dQ_s}{\mathrm dV}.
$$

Compare $C_d$ with depletion capacitance and identify the voltage and frequency regime in which each dominates.

## IV.6 Dynamic prediction

Perform one transient device calculation. Acceptable common-core choices are:

- removal of an injected minority profile after the boundary is switched to equilibrium;
- a forward-to-reverse switching calculation with a declared circuit approximation;
- small-signal response from stored charge; or
- a light-generated pulse collected by the junction.

Report stored charge, boundary current, bulk recombination, and a characteristic switching or collection time. The transient must use the same $D$, $\tau$, geometry, and equilibrium state as the static diode prediction.

## IV.7 Operating limit

Estimate at least one operating limit. Choose one:

- peak field relative to a sourced avalanche or tunneling threshold;
- avalanche voltage from sourced ionization coefficients;
- tunneling onset from a sourced model;
- velocity-saturation limit;
- series-resistance power and self-heating estimate;
- punch-through; or
- a documented reason that the available material data do not support a quantitative limit.

Do not present a silicon empirical coefficient as a prediction for another material.

### Required outputs for Part IV

- material-specific ideal $J$–$V$ curves over temperature;
- analytic-versus-numerical minority profiles and currents;
- ideal and selected nonideal current components;
- extracted ideality factor over bias;
- stored charge and $C_d(V)$;
- comparison of $C_d$ and $C_j$;
- transient profiles and terminal current;
- operating-limit estimate; and
- parameter-sensitivity ranking for $J_s$, forward voltage at fixed current, stored charge, and switching time.

---

# Part V: Integrated prediction, validation, and model boundary

## V.1 Predict one declared device operating point

Select one operating temperature and one practical operating question. Examples include:

- forward voltage required for a specified current density;
- reverse capacitance at a specified voltage;
- Hall voltage and resistivity of either neutral side;
- minority diffusion length and collection fraction;
- stored charge at a specified forward bias;
- recovery time after a declared switching event; or
- maximum reverse voltage before the selected operating-limit model activates.

Report the prediction with an uncertainty interval or sensitivity range. Trace the result backward through the model chain to:

1. LA1 equilibrium inputs;
2. new transport and recombination inputs;
3. geometry and boundary conditions; and
4. numerical approximations.

## V.2 Compare with independent evidence

Compare at least two predicted quantities with independent literature or experimental data. At least one comparison must involve transport or a diode characteristic rather than an equilibrium LA1 quantity.

Acceptable comparisons include:

- mobility, resistivity, or Hall data;
- minority lifetime or diffusion length;
- junction capacitance;
- diode saturation current or forward voltage;
- temperature dependence;
- reverse leakage;
- breakdown or punch-through voltage; or
- switching time.

Do not tune every uncertain parameter to the same curve used for validation. Separate calibration data from validation data.

## V.3 Required validity map

Construct a final validity map with axes chosen from temperature, doping, injection level, bias, field, or device length. Mark regions limited by:

- incomplete ionization;
- degeneracy;
- mobility-model range;
- low-level-injection failure;
- depletion-approximation failure;
- high field or velocity saturation;
- tunneling or avalanche;
- self-heating;
- contact or surface effects; and
- missing material data.

## V.4 Final scientific conclusion

End with a direct answer to:

> To what extent can the LA1 equilibrium material model, augmented by the selected scattering and recombination data, predict the transport and p–n-junction behavior of the declared device?

Your answer must identify:

- the most reliable prediction;
- the least reliable prediction;
- the LA1 parameter with the largest downstream influence;
- the new LA2 parameter with the largest uncertainty;
- the first model assumption that fails as bias or temperature increases; and
- the next experiment or calculation needed to improve the model.

---

## Material-specific pathways

### 4H-SiC or another specified SiC polytype

A defensible SiC project will normally address:

- incomplete donor and acceptor ionization inherited from LA1;
- anisotropic or direction-specific mobility;
- extremely small $n_i$ and its effect on $V_{\mathrm{bi}}$ and ideal $J_s$;
- long minority diffusion scales only when supported by lifetime data;
- high built-in potential and wide depletion regions;
- high critical field as a possible operating advantage; and
- the possibility that measured leakage is controlled by defects or surfaces rather than ideal diffusion.

Log-domain evaluation may be necessary for $n_i$, minority density, and $J_s$. A numerically underflowed zero is not a physical prediction.

### Zinc-blende CdTe

A defensible CdTe project will normally address:

- direct-gap equilibrium statistics inherited from LA1;
- electron and hole mobility data with sample and temperature dependence;
- strong sensitivity of lifetime and junction current to defects and compensation;
- surface recombination and finite thickness where relevant;
- use of CdTe junction or detector literature for validation; and
- distinction between an ideal CdTe homojunction and common heterojunction device structures.

If the available experimental comparison is a CdS/CdTe or other heterojunction, do not compare it directly with a CdTe homojunction without modeling or bounding the band offsets and interface effects.

### Bulk crystalline $\mathrm{Bi_2Te_3}$

This remains an advanced model-validity pathway. A defensible project will normally address:

- anisotropic conductivity and mobility tensors;
- two-carrier Hall response;
- degeneracy and a generalized Einstein relation;
- strong compensation and possible Hall-sign changes;
- narrow-gap reverse leakage and tunneling risk;
- exclusion of topological surface-state transport from the bulk model; and
- the limited validity of a classical abrupt homojunction model.

A p–n-junction calculation is still required, but it may conclude that the textbook depletion and ideal-diode model has a narrow or nonexistent predictive regime. That conclusion must follow from quantitative evidence, not from a general statement that the material is complicated.

---

## Computational organization

The submission must keep inherited data, new data, model code, tests, and generated outputs separate. A recommended structure is

```text
long-assessment-2/
├── README.md
├── requirements.txt
├── data/
│   ├── la1-transfer/
│   │   ├── manifest.yaml
│   │   ├── material-parameters.yaml
│   │   ├── equilibrium-states.csv
│   │   └── source-notes.md
│   ├── transport-parameters.yaml
│   ├── recombination-parameters.yaml
│   ├── device-design.yaml
│   └── validation-data/
├── notebooks/
│   ├── 01_la1_transfer_transport.ipynb
│   ├── 02_recombination_continuity.ipynb
│   ├── 03_pn_electrostatics_cv.ipynb
│   ├── 04_diode_current_switching.ipynb
│   └── 05_integrated_prediction.ipynb
├── src/
│   ├── transport_model.py
│   └── junction_model.py
├── tests/
│   ├── test_transfer.py
│   ├── test_transport_model.py
│   └── test_junction_model.py
├── figures/
└── report.pdf
```

Equivalent organizations are acceptable when the computational argument remains inspectable.

## Minimum reproducibility and validation requirements

The final submission must include automated or reproducible checks for:

1. LA1 transfer-file integrity;
2. charge neutrality of imported equilibrium states;
3. conductivity positivity and correct limiting behavior;
4. one-carrier and two-carrier Hall limits;
5. Einstein-relation units and declared regime;
6. diffusion-length consistency;
7. analytic continuity-equation limits;
8. transient positivity and charge conservation;
9. depletion charge neutrality;
10. Poisson residual and grid convergence;
11. $C_jW=\varepsilon_s$ in the depletion model;
12. zero current at zero junction voltage;
13. analytic versus numerical minority-current agreement;
14. monotonic stored-charge removal after injection is removed; and
15. invariance of physical predictions under a rigid energy-reference shift.

A test that only verifies that code executes is not sufficient physical validation.

## Submission package

Submit:

1. a report of approximately 3,500–5,500 words, excluding equations, captions, tables, and references;
2. the immutable LA1 transfer package;
3. five executed notebooks or an equivalent organized computational record;
4. student-authored reusable Python modules;
5. automated tests and validation tables;
6. all transport, recombination, junction, and validation parameter files;
7. digitized comparison data when used;
8. all plots in vector or high-resolution form;
9. a machine-readable final prediction table; and
10. a short oral defense.

The report must distinguish:

- **inherited state:** quantities established in LA1;
- **new physics:** scattering, recombination, electrostatics, and boundary models;
- **geometry:** device dimensions and junction layout;
- **numerics:** discretization, nonlinear solving, time integration, and convergence;
- **calibration:** data used to determine uncertain parameters;
- **validation:** independent evidence not used for calibration;
- **prediction:** quantities calculated for the declared device; and
- **model boundary:** regimes in which the prediction should not be trusted.

## Oral defense

Prepare a 12-minute defense followed by questions. Include:

- one LA1-to-LA2 model map;
- one parameter-provenance table;
- one numerical validation result;
- one junction diagram;
- one decisive device prediction; and
- one validity map.

Be prepared to answer:

1. Which exact LA1 outputs enter the conductivity calculation?
2. Why can effective mass not determine mobility by itself?
3. What is the difference between drift mobility and Hall mobility?
4. When does $D=\mu k_BT/q$ fail?
5. How did you establish the minority-carrier lifetime?
6. What evidence shows that the continuity solver conserves charge?
7. Why is the built-in potential not simply the band gap divided by $q$?
8. How does incomplete ionization alter depletion width and capacitance?
9. Which LA1 carrier density sets each minority boundary condition?
10. Why does the ideal diode current depend strongly on $n_i$?
11. What causes the first visible departure from ideality in your $J$–$V$ curve?
12. Which parameter dominates uncertainty in the final device prediction?
13. What independent data validate the model?
14. At what bias, field, temperature, or injection level does the model fail first?
15. What result would change if the energy zero in LA1 were shifted by one electron-volt?

## Evaluation approach

The assessment is evaluated holistically as one connected prediction argument. The evaluation will consider:

- integrity of the LA1-to-LA2 transfer;
- correctness and consistency of physical assumptions;
- quality of transport and recombination parameter provenance;
- connection between equilibrium states, transport, and junction behavior;
- dimensional consistency and numerical stability;
- validation against analytic limits and independent data;
- convergence, residual, conservation, and uncertainty evidence;
- distinction between calibrated quantities and genuine predictions;
- interpretation of discrepancies and model limitations;
- clarity of the report, code, plots, and oral defense; and
- the student's ability to defend every link in the prediction chain.

A successful submission does not need to match a commercial device exactly. It must show which behavior follows from the inherited material model, which behavior depends on newly introduced constitutive data, and how uncertainty and model limitations constrain the final prediction.

## Starting references

1. S. M. Sze and M.-K. Lee, *Semiconductor Devices: Physics and Technology*, 3rd ed., Wiley, 2012, Chapters 2–3.
2. S. M. Sze and K. K. Ng, *Physics of Semiconductor Devices*, 3rd ed., Wiley, 2007.
3. A. Schenk, *Advanced Physical Models for Silicon Device Simulation*, Springer, 1998, for model limitations and transport/recombination context.
4. D. K. Schroder, *Semiconductor Material and Device Characterization*, 3rd ed., Wiley, 2006, for mobility, Hall, lifetime, C–V, and diode characterization.
5. The material-specific references and parameter sources established in LA1.
6. Material-specific mobility, lifetime, junction, and high-field literature selected in the LA2 declaration.
