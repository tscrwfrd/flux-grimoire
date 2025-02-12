# FlowForge: Practical Fluid Solver Examples
FlowForge is a comprehensive collection of numerical solvers for computational fluid dynamics (CFD) simulations. This repository serves as a reference implementation for various numerical methods commonly used in fluid dynamics.

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

#### Lax-Wendroff
The [Lax-Wendroff](https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method) The Lax-Wendroff method was introduced by Peter Lax and Burton Wendroff in 1960, and was one of the first second-order numerical schemes developed for solving hyperbolic partial differential equations. So many good references exist, but these are the references I reviewed for this implementation: 
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
- Toro, E. F. (2013). Riemann solvers and numerical methods for fluid dynamics: a practical introduction. Springer Science & Business Media.


#### Roe & Godunov schemes
(TBA)

#### MUSCL
(TBA)

#### WENO
(TBA)


### Running the code

Numerical methods are written in Fortran and plotting routines are implemented in Python's Matplotlib library.  Running the fortran is easiest by using [fpm](https://fpm.fortran-lang.org/) and running:
```bash
fpm run
```
Currently a sucessful run will produce csv file outputs from 1D models in the `data` directory:
- `square_wave.csv` show a square wave being transported with fct
- `dam_break.csv` shows an imaginary dam break using operator splitting
- `sod_shock.csv` shows sod shock problem

Gif animations of model outputs are produced by running the python plotting script. 
Get a virtual environment set-up with [uv](https://docs.astral.sh/uv/) and run:

```bash
uv sync
python python/plotting.py
``` 