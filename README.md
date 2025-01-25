# fluid-forge
This collection implements numerical solvers for fluid dynamics simulations, currently only focusing on the fundamental conservation laws of mass and momentum.

## Fluid Flow Algorithms: Working Examples

This repository contains numerical algorithms I've developed and implemented for various fluid dynamics problems. This repository will continue to grow as I explore new numerical methods and solver implementations.

### Methods

#### Flux-corrected transport
The [Flux Corrected Transport](https://en.wikipedia.org/wiki/Flux-corrected_transport) was a breakthrough in nonlinear finite difference techniques that focused on flux computations.  While the original papers can be challenging to understand, there are other modern references that explain the algorithm clearly, such as:

- Chatterjee, K., & Schunk, R. W. (2020). The development and validation of a ‘flux-corrected transport’based solution methodology for the plasmasphere refilling problem following geomagnetic storms. Earth, Planets and Space, 72, 1-9.

- Kuzmin, D., Löhner, R., & Turek, S. (Eds.). (2012). Flux-corrected transport: principles, algorithms, and applications. Springer Science & Business Media.

**The method is used by propagating a square wave across a structured grid.**

### Running the code

Numerical methods are written in Fortran and plotting routines are implemented in Python's Matplotlib library.  Running the fortran is easiest by using [fpm](https://fpm.fortran-lang.org/) and running:
```bash
fpm run
```
A sucessful run will produce csv file `fct_out.csv` in the root directory.

The easiest way to run the python plotting script is 
getting a virtual environment set-up with [uv](https://docs.astral.sh/uv/) and running:
```bash
uv sync
python python/plot_fct.py
``` 
This produces a gif animation of the data in `fct_out.csv`