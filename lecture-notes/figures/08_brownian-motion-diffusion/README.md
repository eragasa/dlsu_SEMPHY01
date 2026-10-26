# Brownian motion and diffusion notebooks

This directory contains the executed computational sequence supporting the
Brownian-motion and diffusion material. Each notebook separates numerical
evolution, verification, and plotting. Every notebook contains a persisted
inline plot and exports one deterministic PDF figure.

## Student computational labs

The five self-contained labs provide the shortest path from microscopic motion
to macroscopic transport:

| Lab | Computational focus |
|---|---|
| [`lab-01-brownian-to-ficks-law.ipynb`](lab-01-brownian-to-ficks-law.ipynb) | Brownian paths, mean-square displacement, ensemble density, and Fick's first and second laws |
| [`lab-02-particle-current.ipynb`](lab-02-particle-current.ipynb) | Signed crossing counts, particle current, discrete Fick law, local continuity, and charge current |
| [`lab-03-matrix-diffusion-exponential-evolution.ipynb`](lab-03-matrix-diffusion-exponential-evolution.ipynb) | Gradient and Laplacian matrices, Fick flux, matrix-exponential evolution, modal decay, conservation, and explicit-Euler stability |
| [`lab-04-langevin-maxwell-boltzmann-drude.ipynb`](lab-04-langevin-maxwell-boltzmann-drude.ipynb) | Underdamped Langevin dynamics, Maxwell–Boltzmann equilibrium, ballistic-to-diffusive crossover, Einstein relation, and Drude conductivity |
| [`lab-05-effective-mass-from-bands-to-transport.ipynb`](lab-05-effective-mass-from-bands-to-transport.ipynb) | Band-curvature mass tensors, anisotropic Maxwell–Boltzmann velocities, DOS and conductivity masses, multivalley transport, and effective-mass Drude response |

All five notebooks include reproducible random seeds, quantitative assertions,
physical units, interpretation questions, and extension exercises. They are
already executed, but can be rerun from the repository root with:

```bash
uv run --python 3.14 --with jupyter --with numpy --with matplotlib \
  jupyter nbconvert --execute --to notebook --inplace \
  lecture-notes/figures/08_brownian-motion-diffusion/lab-*.ipynb
```

## Supporting computational sequence

| Notebook | Computational focus |
|---|---|
| `01-unbiased-random-walk.ipynb` | Random-walk moments and the mean-square-displacement law |
| `02-wiener-process-quadratic-variation.ipynb` | Wiener refinement, Gaussian increments, and quadratic variation |
| `03-brownian-motion-matrix-covariance.ipynb` | Matrix covariance in one, two, and three dimensions |
| `04-ito-integral-and-formula.ipynb` | Left-endpoint Itô sums and the Itô correction |
| `05-explicit-euler-evolution.ipynb` | Deterministic Euler accuracy and stability |
| `06-euler-maruyama-drift-diffusion.ipynb` | Matrix-valued drift–diffusion simulation |
| `07-sde-to-fokker-planck.ipynb` | Agreement between trajectory ensembles and density evolution |
| `08-master-equation-continuum-limit.ipynb` | Convergence of the discrete master equation to diffusion |
| `09-finite-difference-laplacian.ipynb` | Periodic, reflecting, and absorbing Laplacians |
| `10-diffusive-and-quantum-evolution.ipynb` | Decay versus phase evolution for a shared spatial mode |
| `11-particle-heat-and-carrier-diffusion.ipynb` | Shared operators and distinct transported physical fields |

The corresponding Solid State Physics II copy contains the same supporting
instructional calculations.
