
## The Wave Equation
### Notation
- $|\Psi(t)\rangle$: full quantum state in the Schrödinger picture, depends on time.  
- $|\psi\rangle$: stationary state (eigenvector of $\hat H$), time-independent.  
- $\Psi(x,t) = \langle x|\Psi(t)\rangle$: wavefunction in the $x$-representation.  
- $\psi(x) = \langle x|\psi\rangle$: stationary state wavefunction in the $x$-representation.  
 q     q
##### Time-Dependent Schrödinger Equation
The fundamental postulate is that the state $|\Psi(t)\rangle$ evolves according to
$$
i\hbar \frac{\partial}{\partial t}|\Psi(t)\rangle = \hat H |\Psi(t)\rangle,
$$
where $\hat H$ is the Hamiltonian operator
$$
\hat H = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf{r},t).
$$
$\blacksquare$
### Position Represntation of the TDSE in 1D

Project onto the position basis by multiplying on the left with $\langle x|$:
$$
\langle x|\, i\hbar \frac{\partial}{\partial t}|\Psi(t)\rangle
= \langle x|\hat H|\Psi(t)\rangle.
$$
Define the wavefunction
$$
\Psi(x,t) \equiv \langle x|\Psi(t)\rangle.
$$
Then the TDSE in $x$-representation is
$$
i\hbar \frac{\partial}{\partial t}\Psi(x,t)
= \left(-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x,t)\right)\Psi(x,t).
$$

which gives us
##### Time Dependent Schroedinger equation in 1D
$$
\boxed{
i\hbar \frac{\partial}{\partial t}\Psi(x,t)
= -\frac{\hbar^2}{2m}\frac{\partial^2 \Psi(x,t)}{\partial x^2}
+ V(x,t)\Psi(x,t)}
\tag{TDSE-1D}
$$
$\blacksquare$

### Separability of the TDSE.1D
If the Hamiltonian is time-independent ($V(x,t) = V(x)$), the TDSE is
$$
i\hbar \frac{\partial}{\partial t}\Psi(x,t)
= \left(-\frac{\hbar^2}{2m}\frac{\partial^2}{\partial x^2} + V(x)\right)\Psi(x,t).
$$
The LHS depends on time derivatives while the RHS contains only spatial operators, 
we try a separable form
$$
\Psi(x,t) = \psi(x)\,T(t).
$$

Substituting into the TDSE and dividing by $\psi(x)T(t)$ gives
$$
\frac{i\hbar}{T(t)}\frac{dT(t)}{dt} 
= \frac{1}{\psi(x)}\left(-\frac{\hbar^2}{2m}\frac{d^2}{dx^2}+V(x)\right)\psi(x)
\tag{1}
$$
Since the left-hand side depends only on $t$ and the right-hand side depends only on $x$,  
the equality can hold for all $x$ and $t$ only if both sides are equal to the same constant, we denote this constant as $E$ for energy which will become apparent.

If we set the constant of separation equal to $E$, the spatial equation becomes the 
##### Time Independent Schroedinger Equation in 1D
$$
\boxed{
-\frac{\hbar^2}{2m}\frac{\mathrm{d}^2}{\mathrm{d}x^2}\psi(x) + V(x)\psi(x) = E\psi(x)
}
\tag{TISE.1D}
$$
$\blacksquare$

Equivalently,
$$
\left(-\frac{\hbar^2}{2m}\frac{\mathrm{d}^2}{\mathrm{d}x^2} + V(x)\right)\psi(x) = E\psi(x).
$$
This is an eigenvalue problem where the operator in parentheses is the Hamiltonian $\hat H$ that acts on the stationary state function $\psi(x)$  and $E$ is the eigenvalue.  Since the Hamiltonian $\hat H$ represents the total energy $(\hat{T}+\hat{V}$,  the constant  must be interpreted as the energy eigenvalue of that corresponds with $\psi(x)$.

The LHS of Eqn (1) contain an arbitrary function of $t$, $T(t)$,
$$\begin{align}
   \frac{dT(t)}{dt} = -\frac{iE}{\hbar}\,T(t)
\end{align}$$
Solving this ODE,
##### Time dependence of stationary states
For the TDSE that can be separated by spatial and time components the time evolution of the solutions of the TISE is described by time-dependent phase factor independent of spatial coordinates.
$$
\boxed{T(t) = e^{-iEt/\hbar}}
$$
$\blacksquare$
##### Corollary
Since
$$
|T(t)|^2 = 1,
$$
this phase factor does not affect probabilities or expectation values.  
All measurable physics is contained in $\psi(x)$.  
$\square$

### Evolution of Stationary States
Stationary-state solutions of the TDSE take the form
$$
\Psi(x,t) = \psi(x)\,e^{-iEt/\hbar},
$$
where $\psi(x)$ satisfies the time-independent Schrödinger equation (TISE)
$$
\hat H\psi(x) = E\psi(x),
$$
and $E$ is the corresponding energy eigenvalue.
$\blacksquare$
### Solutions to the TISE in 1D
The spatial function $\psi(x)$ satisfies the eigenvalue problem
$$
\hat{H}\,\psi(x) = E\,\psi(x),
$$
with Hamiltonian operator
$$
\hat H = -\frac{\hbar^2}{2m}\frac{\mathrm{d}^2}{\mathrm{d}x^2} + V(x).
$$
This yields the differential equation
$$
\boxed{
  -\frac{\hbar^2}{2m}\frac{\mathrm{d}^2\psi(x)}{\mathrm{d}x^2}
  + V(x)\,\psi(x) = E\,\psi(x).
}
\tag{TISE.1D}
$$
The solutions take the form  $$
  \{\psi_n(x),\,E_n\}, \quad n=0,1,2,\dots
  $$where each $\psi_n(x)$ is an eigenfunction and $E_n$ is its corresponding eigenvalue, the set $\{\psi_n(x)\}$ forms an orthonormal basis of the Hilbert space $\mathcal{H}$, physical wavefunctions are built from these eigenfunctions according to boundary conditions and normalization.
$$
  \psi(x) = \sum_n c_n\,\psi_n(x)
$$
where $|c_n|^2$ is determined by the occupation statistics.  

For fermions like electrons and nuclei in semiconductors, the relevant distribution is Fermi–Dirac.  For bosonic excitations such as phonons and photons, the occupation is given by Bose–Einstein.  In the non-degenerate (classical) limit $E - \mu \gg k_BT$, both reduce to the Maxwell–Boltzmann distribution.
### General Solution to the TDSE
The full wavefunction is
$$
\Psi(x,t) = \psi(x)\,e^{-iEt/\hbar}.
$$
- The spatial part $\psi(x)$ is determined by boundary conditions and the potential $V(x)$.  
- The probability density is time-independent:  $$
  |\Psi(x,t)|^2 = |\psi(x)|^2.
  $$- Any general state is a superposition of stationary states:
  $$
  \Psi(x,t) = \sum_n c_n\,\psi_n(x)\,e^{-iE_n t/\hbar}.
  $$
## Applications
### Electron in Free Space

#### Infinite Potential Well

#### The Step Potential Function

#### The Potential Barrier

