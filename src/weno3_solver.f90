submodule (fluid_forge) weno3_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  real(real64), parameter :: gamma_g  = 1.4_real64
  real(real64), parameter :: weno_eps = 1.0e-6_real64

  contains

  !> @brief Third-order WENO (Weighted Essentially Non-Oscillatory)
  !> reconstruction with forward-Euler time integration for the 1D
  !> Euler equations.
  !>
  !> Like `roe` and `lax_wendroff`, this is a single-step routine: one
  !> call advances the cell-averaged primitive state (rho, vel, prs)
  !> by one timestep. The spatial discretisation is WENO3; the time
  !> integration is plain forward Euler (first-order in time, like
  !> `roe`). The combination is 1st order in time, formally 3rd order
  !> in space in smooth regions.
  !>
  !> Per step:
  !>   (1) Convert primitives to conserved variables U = (rho, rho*u, E)
  !>       and compute the physical flux F(U) at every cell.
  !>   (2) Global Lax-Friedrichs flux splitting:
  !>           f+_i = (F_i + alpha U_i) / 2,
  !>           f-_i = (F_i - alpha U_i) / 2,
  !>       with alpha = max_i(|u_i| + c_i).
  !>   (3) WENO3 reconstruction of the numerical flux at each face:
  !>           F_face_{i+1/2} = WENO3+(f+_{i-1..i+1})
  !>                          + WENO3-(f-_{i..i+2})
  !>       component-wise on (rho, rho*u, E).
  !>   (4) Forward-Euler conservative update on interior cells:
  !>           U^{n+1}_i = U^n_i - dt/dx * (F_face_{i+1/2} - F_face_{i-1/2})
  !>       converted back to primitives.
  !>   (5) Transmissive (zero-gradient) BCs on the two ghost cells per
  !>       side.
  !>
  !> @note Requires **two ghost cells per side** (4 cells total),
  !>       matching the `roe` / `lax_wendroff` convention. Scenarios can
  !>       reuse nx = 504 directly.
  !> @note Forward Euler with WENO is conditionally stable but does not
  !>       preserve the SSP (Strong Stability Preserving) property of
  !>       the WENO reconstruction. Sharp shocks can produce small
  !>       oscillations; if needed, lower the CFL or wrap this with a
  !>       higher-order RK (the `weno5_rk3_2d` 2D variant uses RK3-SSP
  !>       for that reason).
  !>
  !> Reference for the WENO3 reconstruction (the 3-cell, 2-substencil
  !> classical Jiang-Shu form with optimal weights d_0 = 1/3, d_1 = 2/3
  !> and smoothness indicators beta_k = (Df)^2):
  !>   Jiang, G. S., & Shu, C. W. (1996). Efficient implementation of
  !>     weighted ENO schemes. Journal of Computational Physics, 126(1),
  !>     202-228.   [Section 2 covers WENO3 and WENO5.]
  !>   Shu, C. W. (1998). Essentially non-oscillatory and weighted
  !>     essentially non-oscillatory schemes for hyperbolic conservation
  !>     laws. In Advanced Numerical Approximation of Nonlinear
  !>     Hyperbolic Equations (pp. 325-432). Springer.
  !>     [Lecture-note treatment with the WENO3 weights derivation.]
  !>
  !> @param [inout] rho   Density (nx)
  !> @param [inout] vel   Velocity (nx)
  !> @param [inout] prs   Pressure (nx)
  !> @param [in]    nx    Grid size (including the 4 ghost cells)
  !> @param [in]    dt    Time step
  !> @param [in]    dx    Grid spacing
  module procedure weno3
    integer(int32) :: i, k
    real(real64) :: dtdx, alpha
    real(real64) :: u_loc, c_loc, p_loc, rho_loc
    real(real64) :: u_cons(3, nx), f_phys(3, nx)
    real(real64) :: fp(3, nx), fm(3, nx)
    real(real64) :: f_face(3, nx)

    dtdx = dt / dx

    ! Primitive -> conserved + physical flux + global max wave speed
    alpha = 0.0_real64
    do i = 1, nx
      rho_loc = rho(i)
      u_loc   = vel(i)
      p_loc   = prs(i)
      u_cons(1, i) = rho_loc
      u_cons(2, i) = rho_loc * u_loc
      u_cons(3, i) = p_loc / (gamma_g - 1.0_real64) + 0.5_real64 * rho_loc * u_loc**2
      c_loc = sqrt(gamma_g * p_loc / rho_loc)
      f_phys(1, i) = u_cons(2, i)
      f_phys(2, i) = u_cons(2, i) * u_loc + p_loc
      f_phys(3, i) = u_loc * (u_cons(3, i) + p_loc)
      alpha = max(alpha, abs(u_loc) + c_loc)
    end do

    ! Lax-Friedrichs flux splitting
    do i = 1, nx
      do k = 1, 3
        fp(k, i) = 0.5_real64 * (f_phys(k, i) + alpha * u_cons(k, i))
        fm(k, i) = 0.5_real64 * (f_phys(k, i) - alpha * u_cons(k, i))
      end do
    end do

    ! WENO3 reconstruction at each face i+1/2 (stored at f_face(:, i)).
    ! f+ stencil cells {i-1, i, i+1}; f- stencil cells {i, i+1, i+2}.
    ! Valid range: i in [2, nx-2].
    f_face = 0.0_real64
    do i = 2, nx - 2
      do k = 1, 3
        f_face(k, i) = &
          weno3_recon(fp(k, i-1), fp(k, i),   fp(k, i+1)) + &
          weno3_recon(fm(k, i+2), fm(k, i+1), fm(k, i))
      end do
    end do

    ! Forward-Euler update on interior cells [3, nx-2], then convert
    ! back to primitives.
    do i = 3, nx - 2
      do k = 1, 3
        u_cons(k, i) = u_cons(k, i) - dtdx * (f_face(k, i) - f_face(k, i-1))
      end do
      rho(i) = u_cons(1, i)
      vel(i) = u_cons(2, i) / u_cons(1, i)
      prs(i) = (gamma_g - 1.0_real64) * (u_cons(3, i) - 0.5_real64 * u_cons(2, i)**2 / u_cons(1, i))
      if (prs(i) < 0.0_real64) then
        print *, "Ew...negative pressures being calculated."
      end if
    end do

    ! Transmissive boundary conditions
    rho(1) = rho(3);   rho(2) = rho(3)
    vel(1) = vel(3);   vel(2) = vel(3)
    prs(1) = prs(3);   prs(2) = prs(3)
    rho(nx-1) = rho(nx-2); rho(nx) = rho(nx-2)
    vel(nx-1) = vel(nx-2); vel(nx) = vel(nx-2)
    prs(nx-1) = prs(nx-2); prs(nx) = prs(nx-2)
  end procedure weno3

  !> Classical Jiang-Shu WENO3 reconstruction at face i+1/2 from a
  !> left-biased 3-cell stencil {i-1, i, i+1}.
  !>
  !> The right-biased reconstruction needed for f- at the same face is
  !> obtained by calling with mirrored argument order:
  !>     weno3_recon(g(i+2), g(i+1), g(i))
  !> which exploits the left/right symmetry of the WENO3 formula.
  pure function weno3_recon(fm1, f0, fp1) result(fhat)
    real(real64), intent(in) :: fm1, f0, fp1
    real(real64) :: fhat
    real(real64), parameter :: d0 = 1.0_real64 / 3.0_real64
    real(real64), parameter :: d1 = 2.0_real64 / 3.0_real64
    real(real64) :: q0, q1, beta0, beta1, a0, a1, sum_a

    ! Two 2nd-order candidate reconstructions
    q0 = -0.5_real64 * fm1 + 1.5_real64 * f0
    q1 =  0.5_real64 * f0  + 0.5_real64 * fp1

    ! Smoothness indicators
    beta0 = (f0  - fm1)**2
    beta1 = (fp1 - f0)**2

    ! Non-linear weights
    a0 = d0 / (weno_eps + beta0)**2
    a1 = d1 / (weno_eps + beta1)**2
    sum_a = a0 + a1

    fhat = (a0 * q0 + a1 * q1) / sum_a
  end function weno3_recon

end submodule weno3_impl
