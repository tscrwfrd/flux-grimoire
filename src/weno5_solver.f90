submodule (fluid_forge) weno5_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  real(real64), parameter :: gamma_g  = 1.4_real64
  real(real64), parameter :: weno_eps = 1.0e-6_real64

  contains

  !> @brief Fifth-order WENO (Weighted Essentially Non-Oscillatory)
  !> reconstruction with forward-Euler time integration for the 1D
  !> Euler equations.
  !>
  !> Like `roe` and `lax_wendroff`, one call advances the cell-averaged
  !> primitive state (rho, vel, prs) by one timestep. Spatial scheme is
  !> Jiang-Shu WENO5; the time integration is plain forward Euler. The
  !> combination is 1st order in time, formally 5th order in space in
  !> smooth regions.
  !>
  !> Per step:
  !>   (1) Convert primitives to conserved variables and compute the
  !>       physical flux F(U) at every cell.
  !>   (2) Global Lax-Friedrichs flux splitting:
  !>           f+_i = (F_i + alpha U_i) / 2,
  !>           f-_i = (F_i - alpha U_i) / 2,
  !>       with alpha = max_i(|u_i| + c_i).
  !>   (3) WENO5 reconstruction at each face:
  !>           F_face_{i+1/2} = WENO5+(f+_{i-2..i+2})
  !>                          + WENO5-(f-_{i-1..i+3})
  !>       component-wise on (rho, rho*u, E).
  !>   (4) Forward-Euler conservative update on interior cells,
  !>       converted back to primitives.
  !>   (5) Transmissive BCs on the three ghost cells per side.
  !>
  !> @note Requires **three ghost cells on each side** (6 total), not
  !>       the two used by `roe` / `lax_wendroff` / `weno3`. The WENO5
  !>       5-cell stencil straddles 6 cells away from a cell being
  !>       updated (worst case), and there is no way to lower this
  !>       without changing the reconstruction order. Scenarios use
  !>       nx = 506 for 500 interior cells.
  !> @note Forward Euler with WENO is conditionally stable but does not
  !>       preserve the SSP property; if oscillations appear at strong
  !>       shocks, wrap this with a higher-order RK time scheme..
  !>
  !> References:
  !>   Jiang, G. S., & Shu, C. W. (1996). Efficient implementation of
  !>     weighted ENO schemes. Journal of Computational Physics, 126(1),
  !>     202-228.
  !>   Shu, C. W. (1998). Essentially non-oscillatory and weighted
  !>     essentially non-oscillatory schemes for hyperbolic conservation
  !>     laws. In Advanced Numerical Approximation of Nonlinear
  !>     Hyperbolic Equations (pp. 325-432). Springer.
  !>
  !> @param [inout] rho   Density (nx)
  !> @param [inout] vel   Velocity (nx)
  !> @param [inout] prs   Pressure (nx)
  !> @param [in]    nx    Grid size (including the 6 ghost cells)
  !> @param [in]    dt    Time step
  !> @param [in]    dx    Grid spacing
  module procedure weno5
    integer(int32) :: i, k
    real(real64) :: dtdx, alpha
    real(real64) :: u_loc, c_loc, p_loc, rho_loc
    real(real64) :: u_cons(3, nx), f_phys(3, nx)
    real(real64) :: fp(3, nx), fm(3, nx)
    real(real64) :: f_face(3, nx)

    dtdx = dt / dx

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

    ! WENO5 reconstruction at each face i+1/2 (stored at f_face(:, i)).
    ! f+ stencil cells {i-2, i-1, i, i+1, i+2}; f- stencil cells
    ! {i-1, i, i+1, i+2, i+3}. Valid range: i in [3, nx-3].
    f_face = 0.0_real64
    do i = 3, nx - 3
      do k = 1, 3
        f_face(k, i) = &
          weno5_recon(fp(k, i-2), fp(k, i-1), fp(k, i),   fp(k, i+1), fp(k, i+2)) + &
          weno5_recon(fm(k, i+3), fm(k, i+2), fm(k, i+1), fm(k, i),   fm(k, i-1))
      end do
    end do

    ! Forward-Euler update on interior cells [4, nx-3].
    do i = 4, nx - 3
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

    ! Transmissive BCs on the three ghost cells per side.
    rho(1) = rho(4); rho(2) = rho(4); rho(3) = rho(4)
    vel(1) = vel(4); vel(2) = vel(4); vel(3) = vel(4)
    prs(1) = prs(4); prs(2) = prs(4); prs(3) = prs(4)
    rho(nx-2) = rho(nx-3); rho(nx-1) = rho(nx-3); rho(nx) = rho(nx-3)
    vel(nx-2) = vel(nx-3); vel(nx-1) = vel(nx-3); vel(nx) = vel(nx-3)
    prs(nx-2) = prs(nx-3); prs(nx-1) = prs(nx-3); prs(nx) = prs(nx-3)
  end procedure weno5

  !> Classical Jiang-Shu WENO5 reconstruction at face i+1/2 from a
  !> left-biased 5-cell stencil {i-2, i-1, i, i+1, i+2}.
  pure function weno5_recon(fm2, fm1, f0, fp1, fp2) result(fhat)
    real(real64), intent(in) :: fm2, fm1, f0, fp1, fp2
    real(real64) :: fhat
    real(real64), parameter :: d0 = 0.1_real64, d1 = 0.6_real64, d2 = 0.3_real64
    real(real64), parameter :: r13_12 = 13.0_real64 / 12.0_real64
    real(real64) :: c0, c1, c2, beta0, beta1, beta2, a0, a1, a2, sum_a

    c0 = (2.0_real64 * fm2 - 7.0_real64 * fm1 + 11.0_real64 * f0)  / 6.0_real64
    c1 = (-fm1 + 5.0_real64 * f0 + 2.0_real64 * fp1) / 6.0_real64
    c2 = (2.0_real64 * f0 + 5.0_real64 * fp1 - fp2) / 6.0_real64

    beta0 = r13_12 * (fm2 - 2.0_real64 * fm1 + f0)**2  + &
            0.25_real64 * (fm2 - 4.0_real64 * fm1 + 3.0_real64 * f0)**2
    beta1 = r13_12 * (fm1 - 2.0_real64 * f0 + fp1)**2  + &
            0.25_real64 * (fm1 - fp1)**2
    beta2 = r13_12 * (f0  - 2.0_real64 * fp1 + fp2)**2 + &
            0.25_real64 * (3.0_real64 * f0 - 4.0_real64 * fp1 + fp2)**2

    a0 = d0 / (weno_eps + beta0)**2
    a1 = d1 / (weno_eps + beta1)**2
    a2 = d2 / (weno_eps + beta2)**2
    sum_a = a0 + a1 + a2

    fhat = (a0 * c0 + a1 * c1 + a2 * c2) / sum_a
  end function weno5_recon

end submodule weno5_impl
