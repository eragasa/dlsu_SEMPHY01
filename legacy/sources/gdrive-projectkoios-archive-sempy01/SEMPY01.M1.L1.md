# Module 1: Introduction to the Quantum Theory of Solids
## Lecture 1: Principles of Quantum Mechanics
#### Reading
- [[Neaman3Ed.pdf#page=45|Neaman, Ch 2, Introduction to Quantum Mechanics]]
#### Additional Reading
- [[Hummel4Ed.pdf]]

---
## Principles of Quantum Mechanics

### Quantization of Energy
#### Read:
- [[YF13Ed.pdf#page=1336| Young Freedman, Section 39.5, Continuous Spectra]], sections "Rayleigh and the Ultraviolet Catasophy" and "Planck and the Quantum Hypothesis"

At the start of the 20th century, classical physics failed to explain blackbody radiation.  Planck resolved this by introducing the idea that energy is exchanged in discrete packets, not continuously.  
##### Postulate: Quantization of Energy (Planck, 1900)
The energy of an oscillator of frequency $\nu$ can only take on discrete values:  
$$
E = n h \nu, \quad n = 0,1,2,\dots
$$
- $E$ = allowed energy of the oscillator  
- $n$ = integer quantum number  
- $h = 6.626 \times 10^{-34}\,\mathrm{J \cdot s}$ = Planck’s constant  
- $\nu$ = frequency of oscillation
$\blacksquare$

### Photoelectric Effect (Einstein, 1905)
#### Read:
- [[YF13Ed.pdf#page=1287|Young Freedman, Ch 38.1, Light Absorbed as Photons: The Photoelectric Effect]]
- A. Einstein, “Über einen die Erzeugung und Verwandlung des Lichtes betreffenden heuristischen Gesichtspunkt,” Annalen der Physik, vol. 17, pp. 132–148, 1905. English translation: “On a Heuristic Point of View About the Creation and Conversion of Light,” Wikisource.  Available: https://en.wikisource.org/wiki/Translation:On_a_Heuristic_Point_of_View_about_the_Creation_and_Conversion_of_Light

The energy in a light wave is also contained in discrete packets or bundles.  The particle like packet of energy is a photon, whose energy is given by $E=h\nu$.
$$
	E=h\nu=\frac{hc}{\lambda}
$$
A photon with sufficient energy can knock an electron from the surface of the material.  The minimum amount of energy required to move an electron is called the work function of a material.

in quantum systems all of the energy from one photon is absorbed—if the process is allowed by quantum mechanics—or none at all. Part of the acquired energy is used to liberate the electron from its atomic binding, and the rest contributes to the electron's kinetic energy as a free particle.
#### Work Function $\phi$ 
The energy difference from the Fermi energy and the ionization energy.  In other words, the work function $\phi$ is the energy which is necessary to transport an electron from $E_F$ to infinity.
$$
	\phi=h\nu_0
$$

##### Maximum KE of photoelectron 
$$
	T_{\mathrm{max}}=\frac{1}{2}mv^2=h\nu-h\nu_0
$$
- $T_{\mathrm{max}}$ 
- $h\nu$: the incident photon energy
- $h\nu_0$: the minium energy (or work function) required to remove an electron from the surface 

![[img_maximumKE_Zn.png]]
Source: Wikipedia.  "Photoelectric Effect"

### Wave Particle Duality
- [[YF13Ed.pdf#page=1299|Young Freedman, Section 38.4, Wave Particle Duality, Probability and Uncertainty]]

#### Wave-Particle Duality Principle (de Broglie 1924)
The momentum of a photon $p$ is given by
$$
	p=\frac{h}{\lambda}
$$
Therefore the wavelength of a particle can be expressed as
$$
	\lambda = \frac{h}{p}
$$
where $p$ is the momentum of the particle $\lambda$ is the _de Broglie wavelngth_ of the matter wave.

### Uncertainty Principle (Heisenberg 1927)

##### Uncertainty of momentum and position
It is impossible to simultaneously describe with absolute accuracy the position ($x$) and momentum of a particle ($p$)
$$
	\Delta p \Delta x \geq \hbar
$$
where
$$
	\hbar \equiv\frac{h}{2\pi}=1.054\times10^{-34}\,\mathrm{J\cdot s}
$$
##### Uncertainty of energy and time
It is impossible to simultaneously describe with absolute accuracy the energy of a particle and the instant of time the particle has arrived.
$$
	\Delta E\Delta t \geq \hbar
$$
##### Consequence of the uncertainty principle
We cannot determine the exact position of an electron, instead we have to use _probability density functions_ that allow us to determine the probability that an electron has a particular energy.