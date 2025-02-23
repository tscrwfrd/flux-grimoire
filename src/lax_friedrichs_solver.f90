submodule (fluid_forge) lax_friedrichs_impl

  use iso_fortran_env, only: real64, int32
  implicit none

  contains

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
    module procedure lax_friedrichs
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
      
    end procedure lax_friedrichs

end submodule lax_friedrichs_impl
