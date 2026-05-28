submodule (fluid_forge) fct_2d_solver_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  contains

  !> @brief Zalesak's multidimensional flux-corrected transport for 2D scalar
  !> advection.
  !>
  !> Solves the linear advection equation
  !>
  !>     ∂ρ/∂t + ∂(uρ)/∂x + ∂(vρ)/∂y = 0
  !>
  !> on a uniform 2D grid by blending a low-order monotone flux (donor cell /
  !> upwind) with a high-order flux (Lax-Wendroff). A multidimensional flux
  !> limiter clips the antidiffusive correction so that no new extrema are
  !> introduced in the updated solution.
  !>
  !> Per timestep:
  !>   (1) Donor-cell flux F^L on every x- and y-face.
  !>   (2) Lax-Wendroff flux F^H on every x- and y-face.
  !>   (3) Antidiffusive flux  A = F^H - F^L.
  !>   (4) Transported-diffused solution from the low-order fluxes:
  !>         ρ^td = ρ^n - Δt/Δx (F^L_{i+1/2} - F^L_{i-1/2})
  !>                    - Δt/Δy (F^L_{j+1/2} - F^L_{j-1/2})
  !>   (5) Per-cell extrema ρ^max, ρ^min over a 3x3 stencil of ρ^n ∪ ρ^td.
  !>   (6) Per-cell sums of incoming P+ and outgoing P- antidiffusive flux,
  !>       allowed increase Q+ = ρ^max - ρ^td and decrease Q- = ρ^td - ρ^min,
  !>       and ratios R+- = min(1, Q+-/P+-).
  !>   (7) Per-face limiter coefficient C in [0,1] from the upwind R values.
  !>   (8) Corrected solution:
  !>         ρ^{n+1} = ρ^td - Δt/Δx Δ(CA)_x - Δt/Δy Δ(CA)_y
  !>
  !> @note Requires two ghost cells on each side (4 cells in each dimension).
  !> @note Transmissive (zero-gradient) boundary conditions on rho_adv.
  !> @note Stable in practice for max(|u|Δt/Δx, |v|Δt/Δy) ≲ 0.5.
  !>
  !> Reference:
  !>   Zalesak, S. T. (1979). Fully multidimensional flux-corrected transport
  !>   algorithms for fluids. Journal of Computational Physics, 31(3),
  !>   335-362.
  !>
  !> @param [in]    rho      Scalar field (nx by ny, including ghost cells)
  !> @param [in]    u        x-velocity component (nx by ny)
  !> @param [in]    v        y-velocity component (nx by ny)
  !> @param [in]    nx, ny   Grid dimensions including ghost cells
  !> @param [in]    dtdx     Δt/Δx
  !> @param [in]    dtdy     Δt/Δy
  !> @param [inout] rho_adv  Advected scalar field (nx by ny)
  module procedure fct_2d
    real(real64), parameter :: eps_div = 1.0e-15_real64
    integer(int32) :: i, j
    real(real64) :: ubar, vbar, rho_l, rho_r
    ! Face quantities are stored at the cell-index of the cell to the
    ! left/below the face: flux_x(i, j) lives at face (i+1/2, j) and is
    ! defined for i in [1, nx-1]; flux_y(i, j) lives at (i, j+1/2) and is
    ! defined for j in [1, ny-1]. Unused entries are left at zero.
    real(real64) :: flux_lo_x(nx, ny), flux_lo_y(nx, ny)
    real(real64) :: flux_hi_x(nx, ny), flux_hi_y(nx, ny)
    real(real64) :: a_x(nx, ny), a_y(nx, ny)
    real(real64) :: rho_td(nx, ny), rho_max(nx, ny), rho_min(nx, ny)
    real(real64) :: p_plus(nx, ny), p_minus(nx, ny)
    real(real64) :: q_plus(nx, ny), q_minus(nx, ny)
    real(real64) :: r_plus(nx, ny), r_minus(nx, ny)

    flux_lo_x = 0.0_real64
    flux_lo_y = 0.0_real64
    flux_hi_x = 0.0_real64
    flux_hi_y = 0.0_real64
    a_x = 0.0_real64
    a_y = 0.0_real64
    rho_td = rho
    r_plus = 0.0_real64
    r_minus = 0.0_real64

    ! (1) + (2): x-face fluxes. ubar is the face-averaged velocity.
    do j = 1, ny
      do i = 1, nx - 1
        ubar = 0.5_real64 * (u(i, j) + u(i+1, j))
        rho_l = rho(i, j)
        rho_r = rho(i+1, j)
        flux_lo_x(i, j) = max(ubar, 0.0_real64) * rho_l + &
                          min(ubar, 0.0_real64) * rho_r
        flux_hi_x(i, j) = 0.5_real64 * ubar * (rho_l + rho_r) - &
                          0.5_real64 * dtdx * ubar**2 * (rho_r - rho_l)
        a_x(i, j) = flux_hi_x(i, j) - flux_lo_x(i, j)
      end do
    end do

    ! (1) + (2): y-face fluxes.
    do j = 1, ny - 1
      do i = 1, nx
        vbar = 0.5_real64 * (v(i, j) + v(i, j+1))
        rho_l = rho(i, j)
        rho_r = rho(i, j+1)
        flux_lo_y(i, j) = max(vbar, 0.0_real64) * rho_l + &
                          min(vbar, 0.0_real64) * rho_r
        flux_hi_y(i, j) = 0.5_real64 * vbar * (rho_l + rho_r) - &
                          0.5_real64 * dtdy * vbar**2 * (rho_r - rho_l)
        a_y(i, j) = flux_hi_y(i, j) - flux_lo_y(i, j)
      end do
    end do

    ! (4) Transported-diffused solution from the low-order fluxes.
    do j = 2, ny - 1
      do i = 2, nx - 1
        rho_td(i, j) = rho(i, j) - &
          dtdx * (flux_lo_x(i, j) - flux_lo_x(i-1, j)) - &
          dtdy * (flux_lo_y(i, j) - flux_lo_y(i, j-1))
      end do
    end do

    ! (5) Per-cell extrema over the 3x3 stencil of rho^n ∪ rho^td.
    do j = 2, ny - 1
      do i = 2, nx - 1
        rho_max(i, j) = max( &
          rho(i-1, j-1),    rho(i, j-1),    rho(i+1, j-1),    &
          rho(i-1, j),      rho(i, j),      rho(i+1, j),      &
          rho(i-1, j+1),    rho(i, j+1),    rho(i+1, j+1),    &
          rho_td(i-1, j-1), rho_td(i, j-1), rho_td(i+1, j-1), &
          rho_td(i-1, j),   rho_td(i, j),   rho_td(i+1, j),   &
          rho_td(i-1, j+1), rho_td(i, j+1), rho_td(i+1, j+1))
        rho_min(i, j) = min( &
          rho(i-1, j-1),    rho(i, j-1),    rho(i+1, j-1),    &
          rho(i-1, j),      rho(i, j),      rho(i+1, j),      &
          rho(i-1, j+1),    rho(i, j+1),    rho(i+1, j+1),    &
          rho_td(i-1, j-1), rho_td(i, j-1), rho_td(i+1, j-1), &
          rho_td(i-1, j),   rho_td(i, j),   rho_td(i+1, j),   &
          rho_td(i-1, j+1), rho_td(i, j+1), rho_td(i+1, j+1))
      end do
    end do

    ! (6) Per-cell antidiffusive sums and limiting ratios.
    do j = 2, ny - 1
      do i = 2, nx - 1
        p_plus(i, j) = &
          dtdx * (max(a_x(i-1, j), 0.0_real64) - min(a_x(i, j), 0.0_real64)) + &
          dtdy * (max(a_y(i, j-1), 0.0_real64) - min(a_y(i, j), 0.0_real64))
        p_minus(i, j) = &
          dtdx * (max(a_x(i, j), 0.0_real64) - min(a_x(i-1, j), 0.0_real64)) + &
          dtdy * (max(a_y(i, j), 0.0_real64) - min(a_y(i, j-1), 0.0_real64))

        q_plus(i, j) = rho_max(i, j) - rho_td(i, j)
        q_minus(i, j) = rho_td(i, j) - rho_min(i, j)

        if (p_plus(i, j) > eps_div) then
          r_plus(i, j) = min(1.0_real64, q_plus(i, j) / p_plus(i, j))
        end if
        if (p_minus(i, j) > eps_div) then
          r_minus(i, j) = min(1.0_real64, q_minus(i, j) / p_minus(i, j))
        end if
      end do
    end do

    ! (7) Per-face limiter coefficient. Overwrite a_x/a_y with C*A so the
    ! final update is a simple finite difference.
    do j = 2, ny - 1
      do i = 2, nx - 2
        if (a_x(i, j) >= 0.0_real64) then
          a_x(i, j) = a_x(i, j) * min(r_plus(i+1, j), r_minus(i, j))
        else
          a_x(i, j) = a_x(i, j) * min(r_plus(i, j), r_minus(i+1, j))
        end if
      end do
    end do
    do j = 2, ny - 2
      do i = 2, nx - 1
        if (a_y(i, j) >= 0.0_real64) then
          a_y(i, j) = a_y(i, j) * min(r_plus(i, j+1), r_minus(i, j))
        else
          a_y(i, j) = a_y(i, j) * min(r_plus(i, j), r_minus(i, j+1))
        end if
      end do
    end do

    ! (8) Corrected solution.
    do j = 3, ny - 2
      do i = 3, nx - 2
        rho_adv(i, j) = rho_td(i, j) - &
          dtdx * (a_x(i, j) - a_x(i-1, j)) - &
          dtdy * (a_y(i, j) - a_y(i, j-1))
      end do
    end do

    ! Transmissive boundary conditions on the two ghost rows/columns.
    do j = 1, ny
      rho_adv(1, j)    = rho_adv(3, j)
      rho_adv(2, j)    = rho_adv(3, j)
      rho_adv(nx-1, j) = rho_adv(nx-2, j)
      rho_adv(nx, j)   = rho_adv(nx-2, j)
    end do
    do i = 1, nx
      rho_adv(i, 1)    = rho_adv(i, 3)
      rho_adv(i, 2)    = rho_adv(i, 3)
      rho_adv(i, ny-1) = rho_adv(i, ny-2)
      rho_adv(i, ny)   = rho_adv(i, ny-2)
    end do

  end procedure fct_2d

end submodule fct_2d_solver_impl
