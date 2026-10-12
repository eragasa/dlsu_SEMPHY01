# Introduction to Quantum Mechanics

[[SEMPY01.M1.L1]]

## Assigned Homework
> **Note:** All sketches must be generated using computational tools (Python/Matplotlib, Desmos, MATLAB, etc.).  
> Hand-drawn plots will not be accepted. Annotate axes clearly and include labels for $x$, $V(x,t)$, and parameters.  

Homework for this module should be one .pdf file.   All problems should start on a new page.

If you are using markdown in obsidian, I good way to do page breaks is
``` markdown
---
<div style="page-break-after: always;"></div>
```
but you will not see the html code

## Introduction to Quantum Mechanics

##### Problem 1  
(adapted from Neaman, 3rd Ed., Problem 2.1)  
*Relevance: Standing wave solutions form the basis for quantization in quantum models of solids.*  

The classical wave equation for a two-wire transmission line is given by
$$
\frac{\partial^2 V(x,t)}{\partial x^2} = LC \cdot \frac{\partial^2 V(x,t)}{\partial t^2}.
$$
One possible solution is given by
$$
V(x,t) = \sin(Kx) \cdot \sin(\omega t),
$$
where $K = \frac{n\pi}{a}$ and $\omega = \frac{K}{\sqrt{LC}}$.  
Sketch, on the same graph, the function $V(x,t)$ as a function of $x$ for $0 \leq x \leq a$ and $n=1$ when  
1. $\omega t = 0$  
2. $\omega t = \pi/2$  
3. $\omega t = \pi$  
4. $\omega t = 3\pi/2$  
5. $\omega t = 2\pi$  
##### Problem 2  
(adapted from Neaman, 3rd Ed., Problems 2.2–2.4)  
*Relevance: Traveling wave solutions introduce phase velocity, important for electron and phonon transport in solids.*  

The functions
$$
V(x,t) = \cos\left(\frac{2\pi x}{\lambda} - \omega t\right)
\quad \text{and} \quad
V(x,t) = \cos\left(\frac{2\pi x}{\lambda} + \omega t\right)
$$
are both solutions to the classical wave equation.

1. Sketch, on the same graph, the function $V(x,t)$ as a function of $x$ for $0 \leq x \leq 3\lambda$ when  
   $\omega t = 0, \; 0.25\pi, \; 0.5\pi, \; 0.75\pi, \; \pi$.  

2. Determine the phase velocities of these traveling waves.

#### Principles of Quantum Mechanics
##### Problem 3  
(adapted from Neaman, 3rd Ed., Problem 2.5)  
*Relevance: The photoelectric effect connects work function and threshold wavelength, fundamental to semiconductor surfaces and optoelectronic devices.*  

The work function of a material refers to the minimum energy required to remove an electron from the material.  
Assume that the work function of gold is 4.90 eV and that of cesium is 1.90 eV.  
Calculate the maximum wavelength of light for the photoelectric emission of electrons for gold and cesium.
##### Problem 4  
(adapted from Neaman, 3rd Ed., Problem 2.6)  
Relevance: The de Broglie wavelength of electrons at low energies is fundamental to quantum confinement and tunneling in semiconductors.  

Calculate the de Broglie wavelength, $\lambda = h/p$, for:  
(a) An electron with kinetic energy of (i) 1.0 eV, and (ii) 100 eV.  
(b) A proton with kinetic energy of 1.0 eV.  
(c) A singly ionized tungsten atom with kinetic energy of 1.0 eV.  
(d) A 2000-kg truck traveling at 20 m/s.
##### Problem 5  
(adapted from Neaman, 3rd Ed., Problem 2.7)  
*Relevance: The thermal average energy, momentum, and wavelength of electrons connect directly to carrier statistics in semiconductors.*  

According to classical physics, the average energy of an electron in an electron gas at thermal equilibrium is $3kT/2$.  
Determine, for $T = 300 \, \text{K}$, the average electron energy (in eV), average electron momentum, and the de Broglie wavelength.

##### Problem 6  
(adapted from Neaman, 3rd Ed., Problem 2.8)  
*Relevance: Comparing photon and electron wavelengths is important for understanding optical absorption in semiconductors.*  

An electron and a photon have the same energy.  
At what value of energy (in eV) will the wavelength of the photon be 10 times that of the electron?
##### Problem 7  
(adapted from Neaman, 3rd Ed., Problem 2.9)  
*Relevance: Electron velocity, energy, and wavelength are directly related to carrier transport in semiconductors.*

(a) An electron is moving with a velocity of $2 \times 10^6 \, \text{cm/s}$.  
Determine the electron energy (in eV), momentum, and de Broglie wavelength (in Å).  

(b) The de Broglie wavelength of an electron is 125 Å.  
Determine the electron energy (in eV), momentum, and velocity.
