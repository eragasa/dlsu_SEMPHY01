# Semiconductor Reference Guide

This document introduces the set of semiconductors that will be used in the course project.  For each material, you will find its band gap, common doping systems, applications, typical pitfalls in data collection or modeling, and a difficulty rating (★ out of 5).

## Gallium arsenide (GaAs) ★★☆☆☆
GaAs is a III–V compound with a direct band gap of about 1.42 eV at room temperature. It is heavily used in optoelectronics and high-speed devices. Doping is straightforward: n-type is commonly achieved with Si, Te, or S, while p-type is done with Zn, C, or Be. Applications include LEDs, lasers, and solar cells. The main challenges are that the reported temperature dependence of the band gap varies across sources and that mobility data depends strongly on sample quality. For this reason GaAs is rated at two stars: the literature is abundant, but interpretation requires some care.

## Gallium nitride (GaN) ★★★★☆
GaN has a wide direct band gap of about 3.4 eV and underpins modern blue LEDs and many power electronic devices. n-type doping is usually done with Si or O, while p-type doping with Mg is possible but suffers from incomplete ionization, which complicates analysis. Data can be inconsistent because both wurtzite and cubic phases exist, and the temperature dependence of the band gap is not always reported consistently. GaN earns four stars because it is crucial technologically but messy to model, especially on the p-type side.

## Indium phosphide (InP) ★★★☆☆
InP is another III–V material with a direct band gap of 1.34 eV. Doping is carried out n-type with S, Se, or Si, and p-type with Zn, Cd, or Mg. It finds applications in photodetectors, fiber-optic communication, and high-speed transistors. The main problem is that data on mobility and freeze-out behavior is harder to find, particularly at low temperatures. InP sits in the middle of the scale at three stars: reasonably well covered in the literature, but with gaps that require approximations.

## Zinc oxide (ZnO) ★★★★☆
ZnO is a II–VI semiconductor with a wide band gap around 3.3 eV. It is used in transparent conductors, UV detectors, and piezoelectric devices. n-type doping is achieved with Al, Ga, or In, while p-type doping is unstable and unreliable in practice. The physics is dominated by defects and excitonic effects, and reported data often disagree from one source to another. ZnO is given four stars because of the difficulty in finding reliable parameters and the dominance of defect physics in its behavior.

## Cadmium telluride (CdTe) ★★★☆☆
CdTe has a direct band gap of about 1.5 eV, making it very suitable for thin-film solar cells and radiation detectors. It can be doped n-type with In or Cl and p-type with As, Sb, or P. The major issue is that mobility and carrier lifetimes depend strongly on growth method, so literature values vary widely. The temperature dependence of the band gap is also inconsistently reported. CdTe is rated at three stars: conceptually straightforward, but plagued by sample-dependent data.

## Indium arsenide (InAs) ★★★★★
InAs is a narrow-gap semiconductor with a band gap of only 0.36 eV, which leads to very high intrinsic carrier densities. n-type doping is possible with Si, S, or Te, and p-type with Zn or Be. It is widely used in infrared detectors and ultrafast devices. The main problems are that its high intrinsic concentration makes it hard to model extrinsic doping regimes, and surface states tend to pin the Fermi level. InAs is given the maximum difficulty rating of five stars because of these complications.

## Germanium (Ge) ★★☆☆☆
Ge is an elemental group IV semiconductor with an indirect band gap of 0.66 eV. n-type doping is done with P, As, or Sb, and p-type with B, Ga, or Al. Germanium was historically important in early transistors and is still used in photodetectors and as a high-mobility channel material. The main pitfalls are incomplete ionization at low temperatures and confusion in the literature when Ge data is mixed with Si. For this reason Ge is rated at two stars: generally easy to handle, but with a few traps.

## Silicon carbide (SiC) ★★★★☆
SiC exists in several polytypes (3C, 4H, 6H), with band gaps ranging from 2.3 to 3.3 eV. n-type doping is usually with N or P, and p-type with Al or B. SiC is central to high-power and high-temperature devices. The challenge is that properties vary significantly between polytypes, and data is often reported without specifying which one. This makes interpretation more complex, leading to a rating of four stars.

## If Data Is Missing
Not all of these materials are equally well documented. If data cannot be found, use standard textbook formulas for $n_i(T)$ and $E_g(T)$ with constants, substitute proxy materials only if justified, and clearly state all assumptions. Ranges can be plotted if exact values are unavailable. Always include PDFs of the references used. Missing data will not be penalized, but failure to explain how the gap was handled will be.