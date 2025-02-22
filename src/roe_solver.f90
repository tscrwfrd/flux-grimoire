submodule (fluid_forge) roe_solver_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  contains

    module procedure roe
      real(real64), parameter :: gamma = 1.4
        integer(int32) :: i
        real(real64) :: flux(9), w(6), diff(3)
        real(real64) :: eigvec(3, 3), eigvecinv(3, 3), diag(3,3)
        real(real64) :: dtdx, rhoavg, velavg, enthalpy_avg, enthalpy(2), cs, &
                        denom, alpha, eng
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

    end procedure roe
end submodule roe_solver_impl
