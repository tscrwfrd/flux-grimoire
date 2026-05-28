# Flux-Grimoire: Practical Fluid Solver Examples
Flux-Grimoire is a comprehensive collection of numerical solvers for computational fluid dynamics (CFD) simulations. This repository serves as a reference implementation for various numerical methods commonly used in fluid dynamics.

The primary goals of this repository are to:

- Provide well-documented implementations of fluid dynamics solvers
- Demonstrate different numerical methods for solving fluid flow problems

## Fluid Flow Algorithms: Working Examples

This repository contains numerical algorithms I've developed and implemented for various fluid dynamics problems. This repository will continue to grow as I explore new numerical methods and solver implementations.

### Methods

#### Flux-corrected transport
The [Flux Corrected Transport](https://en.wikipedia.org/wiki/Flux-corrected_transport) was a breakthrough in nonlinear finite difference techniques that focused on flux computations.  While the original papers can be challenging to understand, there are other modern references that explain the algorithm clearly, such as:

- Chatterjee, K., & Schunk, R. W. (2020). The development and validation of a ‘flux-corrected transport’based solution methodology for the plasmasphere refilling problem following geomagnetic storms. Earth, Planets and Space, 72, 1-9.

- Kuzmin, D., Löhner, R., & Turek, S. (Eds.). (2012). Flux-corrected transport: principles, algorithms, and applications. Springer Science & Business Media.

The 1D implementation uses the Boris and Book (1973) SHASTA algorithm with the full flux-limiting procedure described across their series of papers:
- Boris, J.P., & Book, D.L. (1973). Flux-corrected transport. I. SHASTA, a fluid transport algorithm that works. Journal of Computational Physics, 11(1), 38-69.
- Book, D.L., Boris, J.P., & Hain, K. (1975). Flux-corrected transport II: Generalizations of the method. Journal of Computational Physics, 18(3), 248-283.
- Boris, J.P., & Book, D.L. (1976). Flux-corrected transport. III. Minimal-error FCT algorithms. Journal of Computational Physics, 20(4), 397-431.

#### Flux-corrected transport (2D Zalesak)
A two-dimensional extension of FCT using Zalesak's multidimensional flux limiter. This method solves the 2D scalar advection equation ∂ρ/∂t + ∂(uρ)/∂x + ∂(vρ)/∂y = 0 on a uniform grid by blending a low-order monotone flux (donor cell / upwind) with a high-order flux (Lax-Wendroff). The antidiffusive correction is clipped per cell so that no new extrema are introduced — the hallmark of Zalesak's approach.

The implementation follows the 8-step Zalesak algorithm:
1. Compute low-order (donor-cell) fluxes on every x- and y-face
2. Compute high-order (Lax-Wendroff) fluxes on every x- and y-face
3. Form antidiffusive fluxes A = F^H - F^L
4. Obtain the transported-diffused solution ρ^td from the low-order fluxes
5. Compute per-cell extrema ρ^max, ρ^min over a 3×3 stencil
6. Sum incoming and outgoing antidiffusive flux and compute limiting ratios
7. Derive per-face limiter coefficients from the upwind limiting ratios
8. Apply the corrected solution ρ^{n+1} = ρ^td − divergence of the limited antidiffusive flux

The canonical test case is Zalesak's slotted cylinder — a notched disk advected through a rigid-body rotation field. The slot's survival after a full revolution measures the limiter's quality.

Reference:
- Zalesak, S. T. (1979). Fully multidimensional flux-corrected transport algorithms for fluids. Journal of Computational Physics, 31(3), 335-362.

#### Lax-Wendroff
The [Lax-Wendroff](https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method) method was introduced by Peter Lax and Burton Wendroff in 1960, and was one of the first second-order numerical schemes developed for solving hyperbolic partial differential equations. So many good references exist, but these are the references I reviewed for this implementation: 
- Lax, P., & Wendroff, B. (2005). [Systems of conservation laws](https://doi.org/10.1002/cpa.3160130205). In Selected Papers Volume I (pp. 263-283). Springer, New York, NY.
- Toro, E. F. (2013). Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media.
- LeVeque, R. J., & Leveque, R. J. (1992). Numerical methods for conservation laws (Vol. 214). Basel: Birkhäuser.


#### Lax–Friedrichs
The Lax-Friedrichs method is a simple numerical scheme for solving hyperbolic partial differential equations (like wave equations or conservation laws) by replacing the time derivative with a forward difference and the spatial derivative with a centered difference, while averaging neighboring spatial points.

- First-order accurate in time and space
- Introduces numerical diffusion (dissipation)
- Simple to implement but relatively diffusive
- Advantage: Very stable and simple to implement
- Disadvantage: More diffusive than other methods (tends to smooth out solutions)

This method serves as a building block for understanding more sophisticated numerical schemes, though it's often too diffusive for practical applications.

Many good references exist, but these are the references I reviewed for this implementation:
- [Lax-Friedrichs](https://en.wikipedia.org/wiki/Lax%E2%80%93Friedrichs_method) wikipedia page
- Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media.


#### Roe method
The Roe solver is a sophisticated numerical method for solving the Euler equations in computational fluid dynamics. It uses a linearized approximation of the Riemann problem at cell interfaces, incorporating characteristic decomposition to handle wave propagation accurately.

- Second-order accurate in space and time
- Uses eigenvalue decomposition for wave-based solution
- Handles shock waves and contact discontinuities well
- Requires computation of Roe-averaged states
- Implements characteristic decomposition through eigenstructure

Many good references exist, but these are the references I reviewed for this implementation:

- Roe, P. L. (1981). Approximate Riemann solvers, parameter vectors, and difference schemes. Journal of computational physics, 43(2), 357-372.
- Toro, E. F. (2009). Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media.

### Running the code

Numerical methods are written in Fortran and plotting routines are implemented in Python's Matplotlib library.  Running the fortran is easiest by using [fpm](https://fpm.fortran-lang.org/) and running:
```bash
fpm run
```
A successful run will produce five CSV files in the `data` directory:

| File | Solver | Description |
|------|--------|-------------|
| `square_wave.csv` | FCT (1D) | Square wave advection using the Boris-Book flux limiter |
| `dam_break.csv` | FCT (1D) | Dam break with operator splitting for the pressure gradient |
| `sod_shock_lw.csv` | Lax-Wendroff | Sod shock tube (Euler equations) |
| `sod_shock_roe.csv` | Roe | Sod shock tube (Euler equations) |
| `slotted_cylinder.csv` | FCT (2D) | Zalesak's slotted cylinder in a rigid rotation field |

GIF animations of model outputs are produced by running the Python plotting script.
Set up a virtual environment with [uv](https://docs.astral.sh/uv/) (Python 3.11 or 3.12 required) and run:

```bash
fpm run           # generate CSVs first
uv sync           # install matplotlib
python python/plotting.py
``` 