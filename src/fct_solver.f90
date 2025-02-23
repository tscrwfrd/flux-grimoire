submodule (fluid_forge) fct_solver_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  contains

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
  module procedure fct
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

      vel_bnd = 0.5_real64 * [ &
        vel(j-2) + vel(lidx-3), &
        vel(j-1) + vel(j-2), &
        vel(j)   + vel(j-1), &
        vel(j)   + vel(j+1), &
        vel(j+1) + vel(j+2), &
        vel(j+2) + vel(ridx+3) &   
      ]

      rho_bnd = 0.5_real64 * [ &
        rho(j-2) + rho(lidx-3), &
        rho(j-1) + rho(j-2), &
        rho(j)   + rho(j-1), &
        rho(j)   + rho(j+1), &
        rho(j+1) + rho(j+2), &
        rho(j+2) + rho(ridx+3) &  
      ]

      ! flux through cell walls
      flux = rho_bnd * vel_bnd
    
      ! transported quantities
      rho_trns = rho(j-2:j+2) - (dtdx * (flux(2:6) - flux(1:5)))

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
        
  end procedure fct

end submodule fct_solver_impl
