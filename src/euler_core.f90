module euler_core
  !> Shared physics for the 1D Euler solvers.
  !>
  !> Every finite-volume scheme in this library advances the same
  !> conserved state and uses the same physical flux; only the *numerical*
  !> interface flux differs between schemes. Centralising the equation of
  !> state, the conserved vector
  !>
  !>     U = (rho, rho*u, E),     E = p/(gamma - 1) + 1/2 rho u^2
  !>
  !> and the physical flux
  !>
  !>     F(U) = (rho*u, rho*u^2 + p, u(E + p))
  !>
  !> in one place keeps each solver's file focused on the part that maps
  !> to its source paper: the numerical flux F_{i+1/2}.
  use iso_fortran_env, only: real64
  implicit none
  private

  public :: gamma_g, prim_to_cons_flux, flux_from_cons, cons_to_prim

  !> Ratio of specific heats (ideal gas). Single canonical value for all
  !> solvers.
  real(real64), parameter :: gamma_g = 1.4_real64

  contains

    !> Primitive (rho, u, p) at one cell -> conserved U = (rho, rho*u, E)
    !> and physical flux F(U), computed together to share the energy term.
    pure subroutine prim_to_cons_flux(rho, u, p, u_cons, f_phys)
      real(real64), intent(in)  :: rho, u, p
      real(real64), intent(out) :: u_cons(3), f_phys(3)

      u_cons(1) = rho
      u_cons(2) = rho * u
      u_cons(3) = p / (gamma_g - 1.0_real64) + 0.5_real64 * rho * u**2

      f_phys(1) = u_cons(2)
      f_phys(2) = u_cons(2) * u + p
      f_phys(3) = u * (u_cons(3) + p)
    end subroutine prim_to_cons_flux

    !> Physical flux F(U) evaluated directly from a conserved vector.
    !> Needed by predictor-corrector schemes (e.g. Lax-Wendroff) whose
    !> intermediate state is conserved rather than primitive.
    pure function flux_from_cons(u_cons) result(f_phys)
      real(real64), intent(in) :: u_cons(3)
      real(real64) :: f_phys(3)
      real(real64) :: u, p

      u = u_cons(2) / u_cons(1)
      p = (gamma_g - 1.0_real64) * (u_cons(3) - 0.5_real64 * u_cons(2)**2 / u_cons(1))

      f_phys(1) = u_cons(2)
      f_phys(2) = u_cons(2) * u + p
      f_phys(3) = u * (u_cons(3) + p)
    end function flux_from_cons

    !> Conserved U = (rho, rho*u, E) -> primitives (rho, u, p).
    pure subroutine cons_to_prim(u_cons, rho, u, p)
      real(real64), intent(in)  :: u_cons(3)
      real(real64), intent(out) :: rho, u, p

      rho = u_cons(1)
      u   = u_cons(2) / u_cons(1)
      p   = (gamma_g - 1.0_real64) * (u_cons(3) - 0.5_real64 * u_cons(2)**2 / u_cons(1))
    end subroutine cons_to_prim

end module euler_core
