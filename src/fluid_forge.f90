module fluid_forge
  !> A supply of fluid sovers.
  use iso_fortran_env, only: real64, int32
  implicit none
  private

  public :: say_hello, fct, lax_wendroff, lax_friedrichs
contains

  subroutine say_hello
    print *, "Hello, fluid-forge!"
  end subroutine say_hello

  !> @brief Flux Corrected Transport (FCT) algorithm for solving hyperbolic
  !> conservation laws
  !>
  !>
  !> Flux Corrected Transport (FCT) is a numerical method for solving 
  !> hyperbolic partial differential equations that conservatively advects 
  !> quantities across a structured grid while maintaining positivity and 
  !> minimizing numerical diffusion.  This is an explicit implementation 
  !> using the boris and book method.
  !>
  !> @note Requires two ghost cells at each boundary (4 total)
  !> @note Implements transmissive boundary conditions 
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
  !> @param [inout] rho      Density array (length elements)
  !> @param [in]    vel      Velocity array (length elements)
  !> @param [in]    length   Array length (including ghost cells)
  !> @param [in]    dtdx     Time step / grid spacing ratio
  !> @param [inout] rho_adv  Advected density array (length elements)
  pure subroutine fct(rho, vel, length, dtdx, rho_adv)
    implicit none
    integer, intent(in) :: length
    real(real64), intent(in) :: rho(length), vel(length)
    real(real64), intent(inout) :: rho_adv(length)
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
      sgn = sign(1.0_real64, rho_delta(2))
      flux(1) = min(abs(flux(1)), sgn*rho_delta(3), sgn*rho_delta(1))       
      flux(1) = sgn*max(0.0_real64,  flux(1))
      
      sgn = sign(1.0_real64, rho_delta(3))
      flux(2) = min(abs(flux(2)), sgn*rho_delta(2), sgn*rho_delta(4))       
      flux(2) = sgn*max(0.0_real64,  flux(2))

      rho_adv(j) = rho_diffused(3) - flux(2) + flux(1)
      
    end do

    ! boundary conditions
    rho_adv(1) = rho_adv(3);
    rho_adv(2) = rho_adv(3);
    rho_adv(length - 1) = rho_adv(length - 2);
    rho_adv(length) = rho_adv(length - 2);
        
  end subroutine fct


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
  !>   Toro, E. F. (2013). Riemann solvers and numerical methods for fluid 
  !>   dynamics: a practical introduction. Springer Science & Business Media.
  !>
  !> @param [inout] rho    Density array (nx elements)
  !> @param [inout] vel    Velocity array (nx elements)
  !> @param [inout] prs    Pressure array (nx elements)
  !> @param [in]    nx     Number of grid points
  !> @param [in]    dt     Time step size
  !> @param [in]    dx     Grid spacing
  subroutine lax_wendroff(rho, vel, prs, nx, dt, dx)
    integer, intent(in) :: nx
    real(real64), intent(in) :: dt, dx
    real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    real(real64), parameter :: gamma = 1.4
    integer(int32) :: i
    real(real64) :: eng, prsh, flux(9), W(9), U(nx), V(nx), E(nx), &
      uph(3), umh(3), uve(3)
    
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

  end subroutine lax_wendroff

  
  !> @brief Lax-Friedrichs numerical scheme for solving hyperbolic conservation laws
  !>
  !> This Lax-Friedrichs method is first order accurate and has numerical
  !> dissipative and dispersion properties. The method uses a central difference
  !> scheme with averaging of neighboring points:
  !>
  !> U^t+1 = 1/2[ ( U_i+1 +U_i-1) - Δt/Δx(f(U_i+1) - f(U_i-1))  ]
  !>
  !> @note Negative pressures will trigger a warning message
  !> @note Implements transmissive boundary conditions
  !> @note Requires two ghost cells at each boundary (4 total)
  !>
  !> @param [inout] rho    Density array (nx elements)
  !> @param [inout] vel    Velocity array (nx elements) 
  !> @param [inout] prs    Pressure array (nx elements)
  !> @param [in]    nx     Number of grid points
  !> @param [in]    dt     Time step size
  !> @param [in]    dx     Grid spacing
  subroutine lax_friedrichs(rho, vel, prs, nx, dt, dx)
    integer, intent(in) :: nx
    real(real64), intent(in) :: dt, dx
    real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    real(real64), parameter :: gamma = 1.4
    integer(int32) :: i

    real(real64) :: flux(6), w(6), cv(3), dtdx

    dtdx = dt/dx
    
    do i = 3, nx-2
      ! conserved variables rho, momemtum, energy
      w(1) = rho(i-1)
      w(2) = rho(i-1)*vel(i-1)
      w(3) = rho(i-1)*((prs(i-1)/((gamma - 1.0)*rho(i-1))) + 0.5*vel(i-1)**2)
      w(4) = rho(i+1)
      w(5) = rho(i+1)*vel(i+1)
      w(6) = rho(i+1)*((prs(i+1)/((gamma - 1.0)*rho(i+1))) + 0.5*vel(i+1)**2)

      ! j-1
      flux(1) = vel(i-1)*w(1)
      flux(2) = vel(i-1)*w(2) + prs(i-1)
      flux(3) = vel(i-1)*(w(1)*w(3)+prs(i-1))
      ! j+1
      flux(4) = vel(i+1)*w(4)
      flux(5) = vel(i+1)*w(5) + prs(i+1)
      flux(6) = vel(i+1)*(w(4)*w(6)+prs(i+1))

      if(i > 3) then
        rho(i-1) = cv(1)
        vel(i-1) = cv(2)/cv(1)
        prs(i-1) = (gamma - 1.0) * cv(1) * ((cv(3)/cv(1)) - 0.5*vel(i-1)**2)
        if(prs(i-1) < 0.0) then
          print*, "Ew...negative pressures being calculated."
        end if
      end if

      cv(1) = 0.5*(w(1)+w(4) - dtdx*(flux(4)-flux(1)))
      cv(2) = 0.5*(w(2)+w(5) - dtdx*(flux(5)-flux(2)))
      cv(3) = 0.5*(w(3)+w(6) - dtdx*(flux(6)-flux(3)))

      if (i == nx-2) then
        rho(i) = cv(1)
        vel(i) = cv(2)/cv(1)
        prs(i) = (gamma - 1.0) * cv(1) * ((cv(3)/cv(1)) - 0.5*vel(i-1)**2)
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

  end subroutine lax_friedrichs

  subroutine roe(rho, vel, prs, nx, dt, dx)
    integer, intent(in) :: nx
    real(real64), intent(in) :: dt, dx
    real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    real(real64), parameter :: gamma = 1.4
    integer(int32) :: i

    real(real64) :: flux(9), w(6), diff(3)
    real(real64) :: eigvec(3, 3), eigvecinv(3, 3), diag(3,3)
    real(real64) :: dtdx, rhoavg, velavg, enthalpy_avg, enthalpy(2), cs, denom,&
                    alpha, eng

    dtdx = dt/dx

    do i=3, nx-2
      ! Roe averages
      rhoavg = sqrt(rho(i)*rho(i+1))
      denom = sqrt(vel(i)) + sqrt(vel(i+1))
      velavg = sqrt(rho(i))*vel(i) + sqrt(rho(i+1))
      velavg= velavg / denom
      enthalpy(1) = gamma/(gamma-1.0)*prs(i)+0.5*vel(i)**2
      enthalpy(2) = gamma/(gamma-1.0)*prs(i+1)+0.5*vel(i+1)**2
      enthalpy_avg = sqrt(rho(i))*enthalpy(1) + sqrt(rho(i+1))*enthalpy(2)
      enthalpy_avg = enthalpy_avg/denom

      cs = sqrt((gamma-1.0)*(enthalpy_avg - 0.5*velavg**2))

      ! right eigenvectors of jacobian
      eigvec(1,1) = 1.0
      eigvec(2,1) = velavg - cs
      eigvec(3,1) = enthalpy_avg - velavg*cs
      eigvec(1,2) = 1.0
      eigvec(2,2) = velavg
      eigvec(3,2) = 0.5*velavg*velavg
      eigvec(1,3) = 1.0
      eigvec(2,3) = velavg + cs
      eigvec(3,3) = enthalpy_avg + velavg*cs

      !left eigenvector of jacobian
      alpha = (gamma - 1.0)/(2.0*cs)
      eigvecinv(1,1) = 0.5*velavg**2 + (velavg*cs)/(gamma-1.0)             
      eigvecinv(2,1) = (2.0*cs**2)/(gamma - 1.0) - velavg**2
      eigvecinv(3,1) = 0.5*velavg**2 - (velavg*cs)/(gamma-1.0)            
      eigvecinv(1,2) = -velavg - cs/(gamma - 1.0)            
      eigvecinv(2,2) = 2.0*velavg            
      eigvecinv(3,2) = cs/(gamma-1.0) - velavg            
      eigvecinv(1,3) = 1.0            
      eigvecinv(2,3) = -2.0            
      eigvecinv(3,3) = 1.0
      eigvecinv = alpha*eigvecinv

      ! diagonal matrix - eigenvalues
      diag(1,1) = abs(velavg - cs)
      diag(2,1) = 0.0
      diag(3,1) = 0.0
      diag(1,2) = 0.0
      diag(2,2) = abs(velavg)
      diag(3,2) = 0.0
      diag(1,3) = 0.0
      diag(2,3) = 0.0
      diag(3,3) = abs(velavg + cs)
      

      !fluxes
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
      ! j+1
      eng = prs(i+1)/((gamma-1.0)*rho(i+1)) + 0.5*vel(i+1)**2
      flux(7) = rho(i+1)*vel(i+1)
      flux(8) = rho(i+1)*vel(i+1)**2 + prs(i+1)
      flux(9) = vel(i+1)*(rho(i+1)*eng + prs(i+1))

      ! conserved varaibles
      w(1) = rho(i)
      w(2) = rho(i)*vel(i)
      w(3) = rho(i)*((prs(i)/((gamma - 1.0)*rho(i))) + 0.5*vel(i)**2)
      w(4) = rho(i+1)
      w(5) = rho(i+1)*vel(i+1)
      w(6) = rho(i+1)*((prs(i+1)/((gamma - 1.0)*rho(i+1))) + 0.5*vel(i+1)**2)
      diff(1) = w(4) - w(1)
      diff(2) = w(5) - w(2)
      diff(3) = w(6) - w(3)
       
    
    end do

    


  end subroutine roe


  

end module fluid_forge
