module fluid_forge
  !> A supply of fluid sovers.
  use iso_fortran_env, only: real64, int32
  implicit none
  private

  public :: say_hello, fct, lax_wendroff
contains

  subroutine say_hello
    print *, "Hello, fluid-forge!"
  end subroutine say_hello

  !> Flux Corrected Transport (FCT) is a numerical method for solving 
  !> hyperbolic partial differential equations that conservatively advects 
  !> quantities across a structured grid while maintaining positivity and 
  !> minimizing numerical diffusion.  This is an explicit implementation 
  !> using the boris and book method.
  !>
  !> **NOTE**: This implementation assumes arrays are padded with two ghost 
  !> cells at both ends. 
  !>
  !> The algorithm is described across multiple papers:
  !>    Boris, J.P., & Book, D.L. (1973). Flux-corrected transport. I. 
  !>    SHASTA, a fluid transport algorithm that works. Journal of 
  !>    computational physics, 11(1), 38-69.
  !>
  !>    Book, D.L., Boris, J.P., & Hain, K. (1975). Flux-corrected 
  !>    transport II: Generalizations of the method. 
  !>    Journal of Computational Physics, 18(3), 248-283.
  !>
  !>    Boris, J.P., & Book, D.L. (1976). Flux-corrected transport. III. 
  !>    Minimal-error FCT algorithms. 
  !>    Journal of Computational Physics, 20(4), 397-431.
  !>
  !> @param[in] rho - quantity to advect
  !> @param[in] vel - velocities
  !> @param[in] length - array length
  !> @param[in] dtdx - ratio dt/dx
  !> @param[inout] rho_adv - advected rho quantity
  subroutine fct(rho, vel, length, dtdx, rho_adv)
    implicit none
    integer, intent(in) :: length
    real(real64), intent(inout) :: rho(length), vel(length), rho_adv(length)
    real(real64), intent(in) ::  dtdx

    integer(int32) :: j, lidx, ridx
    real(real64), parameter :: r16 = 1.0D0 / 6.0D0, r13 = 1.0D0 / 3.0D0 
    real(real64) :: sgn 
    real(real64) :: rho_bnd(6), vel_bnd(6), flux(6), phi(6), eps(6), nu(6),&
      rho_diffused(5), rho_trns(5), rho_delta(6), mu(2)

    ! cell wall locations index schemes
    ! size 6               size 5            size 4
    ! index 1: -5/2
    ! index 2: -3/2     index 1: -2
    ! index 3: -1/2     index 2: -1       index 1: -3/2
    ! index 4: +1/2     index 3:  0       index 2: -1/2
    ! index 5: +3/2     index 4: +1       index 3: +1/2
    ! index 6: +5/2     index 5: +2       index 4: +3/2
    do j = 3, length - 2

      ! handling boundary cells with padded cells
      if (j == 3) then
        lidx = 4
      else 
        lidx = j 
      end if

      if (j == length-2) then
        ridx = length-3
      else 
        ridx = j 
      end if

      vel_bnd(1:6) = 0.5D0 * (/ vel(j-2) + vel(lidx-3), &
        vel(j-1) + vel(j-2), &
        vel(j)   + vel(j-1), &
        vel(j)   + vel(j+1), &
        vel(j+1) + vel(j+2), &
        vel(j+2) + vel(ridx+3)   /)

      ! Vectorized rho_bnd calculations
      rho_bnd(1:6) = 0.5D0 * (/ rho(j-2) + rho(lidx-3), &
        rho(j-1) + rho(j-2), &
        rho(j)   + rho(j-1), &
        rho(j)   + rho(j+1), &
        rho(j+1) + rho(j+2), &
        rho(j+2) + rho(ridx+3)   /)

      ! flux through cell walls
      flux = rho_bnd * vel_bnd
    
      ! transported quantities
      rho_trns(1:5) = rho(j-2:j+2) - dtdx * (flux(2:6) - flux(1:5))

      ! Δρ - diffusive quantities  
      rho_delta(2:5) = rho(j-1:j+2) - rho(j-2:j+1)
      rho_delta(1) = rho(j-2) - rho(lidx-3)
      rho_delta(6) = rho(ridx+3) - rho(j+2)

      ! Diffusive fluxes
      eps = dtdx * vel_bnd
      ! based on Table 1 in FCT III
      nu = r16 - r16 * eps**2
      flux = nu * rho_delta

      ! Transported-diffused quantities
      rho_diffused(1:5) = rho_trns(1:5) + flux(2:6) - flux(1:5)
     
      ! delta diffused quantities
      rho_delta(1:4) = rho_diffused(2:5) - rho_diffused(1:4)

      ! based on Table 1 in FCT III
      mu(1) = r16 - r16 * eps(3)**2
      mu(2) = r16 - r16 * eps(4)**2

      ! Based on FCT I/II, using transport-diffused quantity
      ! Could also just use the transported quantity (rho_trns)
      ! It's not clear if one is better than the other!
      flux(1) = mu(1) * (rho_diffused(3) - rho_diffused(2))
      flux(2) = mu(2) * (rho_diffused(4) - rho_diffused(3))

      ! anti-diffusive corrections/quantities
      sgn = sign(1.0D0, rho_delta(2))
      flux(1) = min(abs(flux(1)), sgn*rho_delta(3), sgn*rho_delta(1))       
      flux(1) = sgn*max(0.0,  flux(1))
      
      sgn = sign(1.0D0, rho_delta(3))
      flux(2) = min(abs(flux(2)), sgn*rho_delta(2), sgn*rho_delta(4))       
      flux(2) = sgn*max(0.0,  flux(2))

      rho_adv(j) = rho_diffused(3) - flux(2) + flux(1)
      
    end do

    ! boundary conditions
    rho_adv(1) = rho_adv(3);
    rho_adv(2) = rho_adv(3);
    rho_adv(length - 1) = rho_adv(length - 2);
    rho_adv(length) = rho_adv(length - 2);
        
  end subroutine fct

  !> The Lax-Wendroff method is a second-order accurate scheme that can be used 
  !> to solve lienar advection (hyperbolic ∂u/∂t + a∂u/∂x = 0) equation. This 
  !> method suffers from numerical dispersion when dealing with numerical 
  !> discontinuities. The method is captured through two steps:
  !> (Richtmyer)
  !> (1) Half steps 
  !>         Uh_(i+1/2) = 0.5*(U_i+1 + U_i) - 
  !>                             0.5*Δt/Δx(f(U_i+1) - f(U_i))
  !>         Uh_(i-1/2) = 0.5*(U_i + U_i-1) - 
  !>                             0.5*Δt/Δx(f(U_i) - f(U_i-1))
  !> (2) Full step
  !>         U^t+1 = U^t - Δt/Δx(f(Uh_i+1/2) - f(Uh_i-1/2))
  !> 
  !> **NOTE** Two padded cells for the boundaries.
  !> 
  !> 
  !> A great desciption of the method can be found on wikipedia:
  !> https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method
  !> 
  !> and here:
  !>   Toro, E. F. (2013). Riemann solvers and numerical methods for fluid 
  !>   dynamics: a practical introduction. Springer Science & Business Media.
  ! !>   
  ! subroutine lax_wendroff(rho, vel, eng, nx, dt, dx)
  !   integer, intent(in) :: nx
  !   real(real64), intent(in) :: dt, dx
  !   real(real64), intent(inout) :: rho(nx), vel(nx), eng(nx)
  !   real(real64), parameter :: gamma = 1.4
  !   integer(int32) :: i
  !   real(real64) :: prs, flux(9), U(nx), W(nx, 3), V(nx), E(nx), &
  !     uph(3), umh(3)
    
  !   ! conserved quantities
  !   W(:, 1) = rho
  !   W(:, 2) = rho*vel
  !   W(:, 3) = rho*eng

  !   do i = 3, nx - 2
  !     ! Fluxes 
  !     ! j-1
  !     prs = (gamma - 1.0)*rho(i-1)*(eng(i-1) - 0.5*vel(i-1)**2)
  !     flux(1) = rho(i-1)*vel(i-1)
  !     flux(2) = rho(i-1)*vel(i-1)**2 + prs
  !     flux(3) = vel(i-1)*(rho(i-1)*eng(i-1) + prs)

  !     ! j
  !     prs = (gamma - 1.0)*rho(i)*(eng(i) - 0.5*vel(i)**2)
  !     flux(4) = rho(i)*vel(i)
  !     flux(5) = rho(i)*vel(i)**2 + prs
  !     flux(6) = vel(i)*(rho(i)*eng(i) + prs)

  !     !j+1
  !     prs = (gamma - 1.0)*rho(i+1)*(eng(i+1) - 0.5*vel(i+1)**2)
  !     flux(7) = rho(i+1)*vel(i+1)
  !     flux(8) = rho(i+1)*vel(i+1)**2 + prs
  !     flux(9) = vel(i+1)*(rho(i+1)*eng(i+1) + prs)

  !     ! half time step
  !     uph = 0.5*(W(i, :) + W(i+1, :)) -   &
  !           dt/(2.0*dx)*(flux(7:9) - flux(4:6))
  !     umh = 0.5*(W(i, :) + W(i-1, :)) -   &
  !           dt/(2.0*dx)*(flux(4:6) - flux(1:3))

  !     ! flux j + 1/2
  !     uph(2) = uph(2)/uph(1)
  !     uph(3) = uph(3)/uph(1)
  !     prs = (gamma - 1.0)*uph(1)*(uph(3) - 0.5*uph(2)**2)
  !     flux(1) = uph(1)*uph(2)
  !     flux(2) = uph(1)*uph(2)**2 + prs
  !     flux(3) = uph(2)*(uph(1)*uph(3) + prs)
      
  !     ! flux j - 1/2
  !     umh(2) = umh(2)/umh(1)
  !     umh(3) = umh(3)/umh(1)
  !     prs = (gamma - 1.0)*umh(1)*(umh(3) - 0.5*umh(2)**2)
  !     flux(4) = umh(1)*umh(2)
  !     flux(5) = umh(1)*umh(2)**2 + prs
  !     flux(6) = umh(2)*(umh(1)*umh(3) + prs)

  !     ! Full time step
  !     U(i) = W(i, 1) - dt/dx*(flux(1) - flux(4))
  !     V(i) = W(i, 2) - dt/dx*(flux(2) - flux(5))
  !     E(i) = W(i, 3) - dt/dx*(flux(3) - flux(6))

  !   end do    

  !   where(U > 0.0)
  !     rho = U
  !     vel = V/U
  !     eng = E/U
  !   end where
    
  !   ! boundary cells
  !   rho(1) = U(3)
  !   rho(2) = U(3)
  !   vel(1) = V(3)/U(3)
  !   vel(2) = V(3)/U(3)
  !   eng(1) = E(3)/U(3)
  !   eng(2) = E(3)/U(3)
  !   rho(nx) = U(nx-2)
  !   rho(nx-1) = U(nx-2)
  !   vel(nx) = V(nx-2)/U(nx-2)
  !   vel(nx-1) = V(nx-2)/U(nx-2)
  !   eng(nx) = E(nx-2)/U(nx-2)
  !   eng(nx-1) = E(nx-2)/U(nx-2)

  ! end subroutine lax_wendroff

  subroutine lax_wendroff(rho, vel, prs, nx, dt, dx)
    integer, intent(in) :: nx
    real(real64), intent(in) :: dt, dx
    real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    real(real64), parameter :: gamma = 1.4
    integer(int32) :: i
    real(real64) :: eng, prsh, flux(9), U(nx), W(nx, 3), V(nx), E(nx), &
      uph(3), umh(3), Q(3)
    
    ! conserved quantities
    W(:, 1) = rho
    W(:, 2) = rho*vel
    W(:, 3) = rho*((prs/((gamma - 1.0)*rho)) + 0.5*vel**2)

    do i = 3, nx - 2
      ! Fluxes 
      ! j-1
      eng = prs(i-1)/((gamma-1.0)*rho(i-1)) + 0.5*vel(i-1)**2
      flux(1) = rho(i-1)*vel(i-1)
      flux(2) = rho(i-1)*vel(i-1)**2 + prs(i-1)
      flux(3) = vel(i-1)*(rho(i-1)*eng + prs(i-1))

      ! j
      eng = prs(i)/((gamma-1.0)*rho(i)) + 0.5*vel(i)**2
      flux(4) = rho(i)*vel(i)
      flux(5) = rho(i)*vel(i)**2 + prs(i)
      flux(6) = vel(i)*(rho(i)*eng + prs(i))

      !j+1
      eng = prs(i+1)/((gamma-1.0)*rho(i+1)) + 0.5*vel(i+1)**2
      flux(7) = rho(i+1)*vel(i+1)
      flux(8) = rho(i+1)*vel(i+1)**2 + prs(i+1)
      flux(9) = vel(i+1)*(rho(i+1)*eng + prs(i+1))

      ! half time step
      uph = 0.5*(W(i, :) + W(i+1, :)) -   &
            dt/(2.0*dx)*(flux(7:9) - flux(4:6))
      umh = 0.5*(W(i, :) + W(i-1, :)) -   &
            dt/(2.0*dx)*(flux(4:6) - flux(1:3))

      ! flux j + 1/2
      uph(2) = uph(2)/uph(1)
      uph(3) = uph(3)/uph(1)
      prsh = (gamma - 1.0)*uph(1)*(uph(3) - 0.5*uph(2)**2)
      flux(1) = uph(1)*uph(2)
      flux(2) = uph(1)*uph(2)**2 + prsh
      flux(3) = uph(2)*(uph(1)*uph(3) + prsh)
      
      ! flux j - 1/2
      umh(2) = umh(2)/umh(1)
      umh(3) = umh(3)/umh(1)
      prsh = (gamma - 1.0)*umh(1)*(umh(3) - 0.5*umh(2)**2)
      flux(4) = umh(1)*umh(2)
      flux(5) = umh(1)*umh(2)**2 + prsh
      flux(6) = umh(2)*(umh(1)*umh(3) + prsh)

      ! Full time step
      U(i) = W(i, 1) - dt/dx*(flux(1) - flux(4))
      V(i) = W(i, 2) - dt/dx*(flux(2) - flux(5))
      E(i) = W(i, 3) - dt/dx*(flux(3) - flux(6))

    end do    

    where(U > 0.0)
      rho = U
      vel = V/U
      prs = (gamma - 1.0)*rho*((E/U) - 0.5*vel**2)
    end where
    
    ! boundary cells
    rho(1) = rho(3)
    rho(2) = rho(3)
    vel(1) = vel(3)
    vel(2) = vel(3)
    prs(1) = (gamma - 1.0)*rho(3)*((E(3)/U(3)) - 0.5*vel(3)**2)
    prs(2) = prs(1)

    rho(nx) = rho(nx-2)
    rho(nx-1) = rho(nx-2)
    vel(nx) = vel(nx-2)
    vel(nx-1) = vel(nx)
    prs(nx) = (gamma - 1.0)*rho(nx-2)*((E(nx-2)/U(nx-2)) - 0.5*(V(nx-2)/U(nx-2))**2)
    prs(nx-1) = prs(nx)

  end subroutine lax_wendroff

end module fluid_forge
