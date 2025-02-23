submodule (fluid_forge) lax_wendroff_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  contains
    !> @brief Lax-Wendroff scheme for solving hyperbolic conservation laws
    !>
    !> The Lax-Wendroff method is a second-order accurate scheme that can be used 
    !> to solve lienar advection (hyperbolic ∂u/∂t + a∂u/∂x = 0) equation. This 
    !> method suffers from numerical dispersion when dealing with numerical 
    !> discontinuities. The method is captured through two steps:
    !> (Richtmyer variant)
    !> (1) Half steps 
    !>         Uh_(i+1/2) = 0.5*(U_i+1 + U_i) - 
    !>                             0.5*Δt/Δx(f(U_i+1) - f(U_i))
    !>         Uh_(i-1/2) = 0.5*(U_i + U_i-1) - 
    !>                             0.5*Δt/Δx(f(U_i) - f(U_i-1))
    !> (2) Full step
    !>         U^t+1 = U^t - Δt/Δx(f(Uh_i+1/2) - f(Uh_i-1/2))
    !> 
    !> @note Requires two ghost cells at each boundary (4 total)
    !> @note Implements transmissive boundary conditions
    !> @note Checks for and warns about negative pressures
    !> @note Method exhibits numerical dispersion near discontinuities
    !> 
    !> A great desciption of the method can be found on wikipedia:
    !> https://en.wikipedia.org/wiki/Lax%E2%80%93Wendroff_method
    !> 
    !> and here:
    !>   Toro, E. F. (2009). Riemann solvers and numerical methods for fluid 
    !>   dynamics: a practical introduction. Springer Science & Business Media.
    !>
    !> @param [inout] rho    Density array (nx elements)
    !> @param [inout] vel    Velocity array (nx elements)
    !> @param [inout] prs    Pressure array (nx elements)
    !> @param [in]    nx     Number of grid points
    !> @param [in]    dt     Time step size
    !> @param [in]    dx     Grid spacing
    module procedure lax_wendroff
      real(real64), parameter :: gamma = 1.4
      integer(int32) :: i
      real(real64) :: eng, prsh, flux(9), W(9), uph(3), umh(3), uve(3)

      do i = 3, nx - 2
      ! conserved quantities
      W(1) = rho(i-1)
      W(2) = rho(i-1)*vel(i-1)
      W(3) = rho(i-1)*((prs(i-1)/((gamma - 1.0)*rho(i-1))) + 0.5*vel(i-1)**2)
      W(4) = rho(i)
      W(5) = rho(i)*vel(i)
      W(6) = rho(i)*((prs(i)/((gamma - 1.0)*rho(i))) + 0.5*vel(i)**2)
      W(7) = rho(i+1)
      W(8) = rho(i+1)*vel(i+1)
      W(9) = rho(i+1)*((prs(i+1)/((gamma - 1.0)*rho(i+1))) + 0.5*vel(i+1)**2)


      ! Fluxes 
      ! j-1
      eng = prs(i-1)/((gamma-1.0)*rho(i-1)) + 0.5*vel(i-1)**2
      flux(1) = rho(i-1)*vel(i-1)
      flux(2) = rho(i-1)*vel(i-1)**2 + prs(i-1)
      flux(3) = vel(i-1)*(rho(i-1)*eng + prs(i-1))

      ! now update the solution for location i-1
      if (i > 3) then
        rho(i-1) = uve(1)
        vel(i-1) = uve(2)/uve(1)
        prs(i-1) = (gamma - 1.0)*uve(1)*((uve(3)/uve(1)) - 0.5*vel(i-1)**2)
        if(prs(i-1) < 0.0) then
          print*, "Ew...negative pressures being calculated."
        end if
      end if 

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
      uph = 0.5*(W(4:6) + W(7:9)) -   &
            dt/(2.0*dx)*(flux(7:9) - flux(4:6))
      umh = 0.5*(W(4:6) + W(1:3)) -   &
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
      uve = W(4:6) - dt/dx*(flux(1:3) - flux(4:6))

      ! update since this is last location to compute
      if (i==nx-2) then
        rho(i) = uve(1)
        vel(i) = uve(2)/uve(1)
        prs(i) = (gamma - 1.0) * uve(1) * ((uve(3)/uve(1)) - 0.5*vel(i)**2)
        if(prs(i) < 0.0) then
          print*, "Ew...negative pressures being calculated."
        end if
      end if

    end do    
    
    ! boundary cells
    rho(1) = rho(3)
    rho(2) = rho(3)
    vel(1) = vel(3)
    vel(2) = vel(3)
    prs(1) = prs(3)
    prs(2) = prs(3)

    rho(nx) = rho(nx-2)
    rho(nx-1) = rho(nx)
    vel(nx) = vel(nx-2)
    vel(nx-1) = vel(nx)
    prs(nx) = prs(nx-2)
    prs(nx-1) = prs(nx)
    
    end procedure lax_wendroff

end submodule lax_wendroff_impl
