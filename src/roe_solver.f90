submodule (fluid_forge) roe_solver_impl
  use iso_fortran_env, only: real64, int32
  implicit none

  contains

    module procedure roe
      real(real64), parameter :: gamma = 1.4
      integer(int32) :: i
      real(real64) :: flux(9), w(9), diff(3), qnty(3), update(3), roe_flux(3, 2)
      real(real64) :: eigvec(3, 3), eigvecinv(3, 3), diag(3,3), roe_jac(3,3)
      real(real64) :: dtdx, eng
      dtdx = dt/dx

      do i=3, nx-2

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

        ! conserved variables
        w(1) = rho(i-1)
        w(2) = rho(i-1)*vel(i-1)
        w(3) = rho(i-1)*((prs(i-1)/((gamma - 1.0)*rho(i-1))) + 0.5*vel(i-1)**2)
        w(4) = rho(i)
        w(5) = rho(i)*vel(i)
        w(6) = rho(i)*((prs(i)/((gamma - 1.0)*rho(i))) + 0.5*vel(i)**2)
        w(7) = rho(i+1)
        w(8) = rho(i+1)*vel(i+1)
        w(9) = rho(i+1)*((prs(i+1)/((gamma - 1.0)*rho(i+1))) + 0.5*vel(i+1)**2)
        
        call jac_decomp(rho, vel, prs, i-1, eigvec, eigvecinv, diag)
        roe_jac = matmul(eigvec, diag)
        roe_jac = matmul(roe_jac, eigvecinv)
        diff(1) = w(7) - w(4)
        diff(2) = w(8) - w(5)
        diff(3) = w(9) - w(6)
        qnty = matmul(roe_jac, diff)
        ! Flux_{j-1/2}
        roe_flux(:, 1) = 0.5_real64*(flux(1:3) + flux(4:6)) - 0.5_real64*(qnty)

        call jac_decomp(rho, vel, prs, i, eigvec, eigvecinv, diag)
        roe_jac = matmul(eigvec, diag)
        roe_jac = matmul(roe_jac, eigvecinv)
        diff(1) = w(4) - w(1)
        diff(2) = w(5) - w(2)
        diff(3) = w(6) - w(3)
        qnty = matmul(roe_jac, diff)
        ! Flux_{j+1/2}
        roe_flux(:, 2) = 0.5_real64*(flux(4:6) + flux(7:9)) - 0.5_real64*(qnty)
        
        ! scheme update: U^t+1 = U^t - Δt/Δx(f(U_i+1) - f(U_i-1)
        w(4:6) = w(4:6) - dtdx*(roe_flux(:, 2) - roe_flux(:, 1))

        ! if (i == 5) then
        !   print *, rho(i), qnty
        ! end if

        ! update last
        if (i > 3) then 
          rho(i-1) = update(1)
          vel(i-1) = update(2)/update(1)
          prs(i-1) = (gamma - 1.0_real64)*update(1)*((update(3)/update(1)) - 0.5_real64*vel(i-1)**2)
          if(prs(i-1) < 0.0) then
            print*, "Ew...negative pressures being calculated."
          end if
        end if

        update = w(4:6)

        if (i == nx-2) then 
          rho(i) = update(1)
          vel(i) = update(2)/update(1)
          prs(i) = (gamma - 1.0_real64)*update(1)*((update(3)/update(1)) - 0.5_real64*vel(i)**2)
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

    end procedure roe

    subroutine jac_decomp(rho, vel, prs, idx, eigvec, eigvecinv, diag)
      real(real64), intent(in) :: rho(:), vel(:), prs(:)
      integer(int32), intent(in) :: idx
      real(real64), intent(out) :: eigvec(3,3), eigvecinv(3,3), diag(3,3)
      real(real64), parameter :: gamma = 1.4
      real(real64) :: rhoavg, velavg, enthalpy_avg, enthalpy(2), cs, &
                        denom, alpha
      
      ! Roe averages
      rhoavg = sqrt(rho(idx)*rho(idx+1))
      denom = sqrt(vel(idx)) + sqrt(vel(idx+1))
      velavg = sqrt(rho(idx))*vel(idx) + sqrt(rho(idx+1))
      velavg= velavg / denom
      enthalpy(1) = gamma/(gamma-1.0)*prs(idx)+0.5*vel(idx)**2
      enthalpy(2) = gamma/(gamma-1.0)*prs(idx+1)+0.5*vel(idx+1)**2
      enthalpy_avg = sqrt(rho(idx))*enthalpy(1) + sqrt(rho(idx+1))*enthalpy(2)
      enthalpy_avg = enthalpy_avg/denom

      cs = sqrt((gamma-1.0)*(enthalpy_avg - 0.5*velavg**2))
      if (idx == 5) then
        print *, rho(idx)
      end if

      ! right eigenvectors of jacobian
      eigvec(1,1) = 1.0
      eigvec(2,1) = velavg - cs
      eigvec(3,1) = enthalpy_avg - velavg*cs
      eigvec(1,2) = 1.0
      eigvec(2,2) = velavg
      eigvec(3,2) = 0.5*velavg**2
      eigvec(1,3) = 1.0
      eigvec(2,3) = velavg + cs
      eigvec(3,3) = enthalpy_avg + velavg*cs

      !left eigenvector of jacobian
      alpha = (gamma - 1.0)/(2.0*cs**2)
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


    end subroutine jac_decomp

end submodule roe_solver_impl
