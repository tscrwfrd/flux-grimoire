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

#### WENO3 (third-order)
The Weighted Essentially Non-Oscillatory (WENO) method is a high-order accurate scheme that avoids Gibbs phenomena at discontinuities by adaptively weighting several candidate stencils. This 1D implementation uses the classical Jiang-Shu WENO3 reconstruction with global Lax-Friedrichs flux splitting and forward-Euler time integration. The combination is 1st-order in time, formally 3rd-order in space in smooth regions.

Per step:
1. Convert primitives to conserved variables U = (ρ, ρu, E) and compute the physical flux F(U) at every cell
2. Global Lax-Friedrichs flux splitting: f⁺ = (F + αU)/2, f⁻ = (F − αU)/2, with α = max(|u| + c)
3. WENO3 reconstruction of the numerical flux at each face using two 2nd-order candidate substencils (optimal weights d₀ = 1/3, d₁ = 2/3) with Jiang-Shu smoothness indicators βₖ = (Δf)²
4. Forward-Euler conservative update on interior cells, converted back to primitives
5. Transmissive (zero-gradient) BCs on the two ghost cells per side

Requires **two ghost cells per side** (4 total, matching the `roe` / `lax_wendroff` convention), so scenarios use nx = 504 for 500 interior cells. Forward Euler with WENO is conditionally stable but does not preserve the SSP property of the spatial reconstruction; sharp shocks can produce small oscillations.

References:
- Jiang, G. S., & Shu, C. W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228.
- Shu, C. W. (1998). Essentially non-oscillatory and weighted essentially non-oscillatory schemes for hyperbolic conservation laws. In Advanced Numerical Approximation of Nonlinear Hyperbolic Equations (pp. 325-432). Springer.

#### WENO5 (fifth-order)
A fifth-order extension of the same approach. WENO5 uses a 5-cell stencil with three 3rd-order candidate substencils and the classical Jiang-Shu weights (d₀ = 0.1, d₁ = 0.6, d₂ = 0.3) for the reconstruction at each face. The time integration is forward Euler (1st-order in time, formally 5th-order in space in smooth regions).

The per-step structure mirrors WENO3, except:
- The WENO5 stencil spans five cells, requiring **three ghost cells on each side** (6 total) — one more than the other solvers. Scenarios use nx = 506 for 500 interior cells.
- The reconstruction at face i+1/2 uses f⁺ on cells {i−2, …, i+2} and f⁻ on cells {i−1, …, i+3}.

Forward Euler with WENO is conditionally stable but does not preserve the SSP property; if oscillations appear at strong shocks, the spatial scheme should be wrapped with RK3-SSP.

References:
- Jiang, G. S., & Shu, C. W. (1996). Efficient implementation of weighted ENO schemes. Journal of Computational Physics, 126(1), 202-228.
- Shu, C. W. (1998). Essentially non-oscillatory and weighted essentially non-oscillatory schemes for hyperbolic conservation laws. In Advanced Numerical Approximation of Nonlinear Hyperbolic Equations (pp. 325-432). Springer.

### Code architecture

The library is built around a parent module, `fluid_forge`, that declares the public interface for every solver. Each scheme lives in its own **submodule** so that a single source file maps cleanly to the part of the algorithm that corresponds to its source paper:

| File | Contains |
|------|----------|
| `src/fluid_forge.f90` | Parent module: public interface declarations (`fct`, `fct_2d`, `lax_wendroff`, `lax_friedrichs`, `roe`, `weno3`, `weno5`) |
| `src/euler_core.f90` | **Shared physics for the 1D Euler solvers** (see below) |
| `src/lax_wendroff_solver.f90` | Lax-Wendroff (Richtmyer two-step) submodule |
| `src/lax_friedrichs_solver.f90` | Lax-Friedrichs submodule |
| `src/roe_solver.f90` | Roe approximate Riemann solver + Roe-average eigendecomposition |
| `src/weno3_solver.f90` | WENO3 reconstruction submodule |
| `src/weno5_solver.f90` | WENO5 reconstruction submodule |
| `src/fct_solver.f90`, `src/fct_2d_solver.f90` | 1D and 2D flux-corrected transport |
| `src/fluid_1d_models.f90`, `src/fluid_2d_models.f90` | Scenario drivers that set up each test case and write the CSV output |

#### Shared Euler physics (`euler_core`)

Every finite-volume scheme for the Euler equations advances the same conserved state and uses the same physical flux; only the *numerical* interface flux F_{i+1/2} differs between schemes. `euler_core` centralises the parts that are common to all of them:

- `gamma_g` — the single canonical ratio of specific heats (γ = 1.4) used by every Euler solver.
- `prim_to_cons_flux(rho, u, p, u_cons, f_phys)` — primitives → conserved vector `U = (ρ, ρu, E)` **and** physical flux `F(U) = (ρu, ρu² + p, u(E + p))`, computed together to share the energy term.
- `flux_from_cons(u_cons)` — physical flux evaluated directly from a conserved vector, for predictor-corrector schemes (e.g. Lax-Wendroff) whose intermediate state is conserved rather than primitive.
- `cons_to_prim(u_cons, rho, u, p)` — conserved vector → primitives.

The `lax_wendroff`, `roe`, `weno3`, and `weno5` submodules all `use euler_core` for these operations, so the equation of state lives in exactly one place. Each solver's file is left focused on the numerical flux that maps to its source paper — the Roe-average eigenstructure, the WENO reconstruction weights, and so on.

### Running the code

Numerical methods are written in Fortran and plotting routines are implemented in Python's Matplotlib library.  Running the fortran is easiest by using [fpm](https://fpm.fortran-lang.org/) and running:
```bash
fpm run
```
A successful run will produce seven CSV files in the `data` directory:

| File | Executable | Solver | Description |
|------|------------|--------|-------------|
| `square_wave.csv` | `square-wave` | FCT (1D) | Square wave advection using the Boris-Book flux limiter |
| `dam_break.csv` | `dam-break` | FCT (1D) | Dam break with operator splitting for the pressure gradient |
| `sod_shock_lw.csv` | `sod-shock-lw` | Lax-Wendroff | Sod shock tube (Euler equations) |
| `sod_shock_roe.csv` | `sod-shock-roe` | Roe | Sod shock tube (Euler equations) |
| `slotted_cylinder.csv` | `slotted-cylinder` | FCT (2D) | Zalesak's slotted cylinder in a rigid rotation field |
| `sod_shock_weno3.csv` | `sod-shock-weno3` | WENO3 | Sod shock tube (WENO3 reconstruction, forward Euler) |
| `sod_shock_weno5.csv` | `sod-shock-weno5` | WENO5 | Sod shock tube (WENO5 reconstruction, forward Euler) |

#### Running an individual model

Each model also builds as its own executable, so you can run any single one without running the rest. Pass the executable name (from the table above) to `fpm run`:

```bash
fpm run sod-shock-weno5      # run just the WENO5 Sod shock tube
fpm run slotted-cylinder     # run just the 2D slotted cylinder
```

To see every available target:

```bash
fpm run --list
```

GIF animations of model outputs are produced by running the Python plotting script.
Set up a virtual environment with [uv](https://docs.astral.sh/uv/) (Python 3.11 or 3.12 required) and run:

```bash
fpm run           # generate CSVs first
uv sync           # install matplotlib
python python/plotting.py
``` 