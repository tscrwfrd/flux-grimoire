module fluid_forge
  !> A supply of fluid sovers.
  use iso_fortran_env, only: real64, int32
  implicit none
  private

  public :: say_hello, fct, lax_wendroff, lax_friedrichs, roe

  interface

    module subroutine fct(rho, vel, length, dtdx, rho_adv)
      integer, intent(in) :: length
      real(real64), intent(in) :: rho(length), vel(length)
      real(real64), intent(inout) :: rho_adv(length)
      real(real64), intent(in) ::  dtdx
    end subroutine fct

    module subroutine lax_wendroff(rho, vel, prs, nx, dt, dx)
      integer, intent(in) :: nx
      real(real64), intent(in) :: dt, dx
      real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    end subroutine lax_wendroff

    module subroutine lax_friedrichs(rho, vel, prs, nx, dt, dx)
      integer, intent(in) :: nx
      real(real64), intent(in) :: dt, dx
      real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    end subroutine lax_friedrichs

    module subroutine roe(rho, vel, prs, nx, dt, dx)
      integer, intent(in) :: nx
      real(real64), intent(in) :: dt, dx
      real(real64), intent(inout) :: rho(nx), vel(nx), prs(nx)
    end subroutine roe
      
  end interface

  contains

    subroutine say_hello
      print *, "Hello, fluid-forge!"
    end subroutine say_hello

end module fluid_forge
