program check
  !> Unit tests for euler_core, the shared 1D-Euler physics used by the
  !> lax_wendroff, roe, weno3 and weno5 solvers. These exercise the
  !> equation of state and flux routines directly so a regression in the
  !> shared core is caught here rather than as a subtle drift in every
  !> solver's CSV output.
  use iso_fortran_env, only: real64, error_unit
  use euler_core, only: gamma_g, prim_to_cons_flux, flux_from_cons, cons_to_prim
  implicit none

  integer :: nfail = 0
  real(real64), parameter :: tol = 1.0e-9_real64

  call test_gamma()
  call test_known_values()
  call test_zero_velocity()
  call test_prim_cons_roundtrip()
  call test_flux_consistency()

  if (nfail == 0) then
    print *, "All euler_core tests passed."
  else
    write(error_unit, '(A,I0,A)') "FAILED: ", nfail, " euler_core check(s) failed."
    error stop 1
  end if

contains

  !> Compare two scalars with a mixed absolute/relative tolerance and
  !> record a failure (without aborting) so every check in a test runs.
  subroutine assert_close(got, want, name)
    real(real64), intent(in) :: got, want
    character(*), intent(in) :: name
    if (abs(got - want) > tol * max(1.0_real64, abs(want))) then
      write(error_unit, '(A,A,A,ES23.15,A,ES23.15)') &
        "  FAIL ", name, ": got ", got, ", want ", want
      nfail = nfail + 1
    end if
  end subroutine assert_close

  !> The single canonical ratio of specific heats. The hand-computed
  !> expected values below assume gamma = 1.4; this guards that premise.
  subroutine test_gamma()
    call assert_close(gamma_g, 1.4_real64, "gamma_g")
  end subroutine test_gamma

  !> Hand-computed conserved vector and physical flux for a known state.
  !> rho=2, u=3, p=5, gamma=1.4:
  !>   E       = p/(gamma-1) + 0.5*rho*u^2 = 12.5 + 9 = 21.5
  !>   U       = (2, 6, 21.5)
  !>   F(U)    = (rho*u, rho*u^2 + p, u*(E + p)) = (6, 23, 79.5)
  subroutine test_known_values()
    real(real64) :: u_cons(3), f_phys(3)
    call prim_to_cons_flux(2.0_real64, 3.0_real64, 5.0_real64, u_cons, f_phys)
    call assert_close(u_cons(1),  2.0_real64,  "known U(1)")
    call assert_close(u_cons(2),  6.0_real64,  "known U(2)")
    call assert_close(u_cons(3), 21.5_real64,  "known U(3)")
    call assert_close(f_phys(1),  6.0_real64,  "known F(1)")
    call assert_close(f_phys(2), 23.0_real64,  "known F(2)")
    call assert_close(f_phys(3), 79.5_real64,  "known F(3)")
  end subroutine test_known_values

  !> At rest (u=0) the momentum flux is pure pressure and the mass and
  !> energy fluxes vanish; internal energy is p/(gamma-1).
  subroutine test_zero_velocity()
    real(real64) :: u_cons(3), f_phys(3)
    call prim_to_cons_flux(1.0_real64, 0.0_real64, 1.0_real64, u_cons, f_phys)
    call assert_close(u_cons(2), 0.0_real64, "rest U(2)")
    call assert_close(u_cons(3), 1.0_real64 / (gamma_g - 1.0_real64), "rest U(3)")
    call assert_close(f_phys(1), 0.0_real64, "rest F(1)")
    call assert_close(f_phys(2), 1.0_real64, "rest F(2)")
    call assert_close(f_phys(3), 0.0_real64, "rest F(3)")
  end subroutine test_zero_velocity

  !> cons_to_prim must invert prim_to_cons_flux: primitives -> conserved
  !> -> primitives recovers the original state. Uses a realistic
  !> (non-trivial) state to exercise the energy/velocity round trip.
  subroutine test_prim_cons_roundtrip()
    real(real64), parameter :: rho0 = 1.225_real64, u0 = 37.0_real64, p0 = 101325.0_real64
    real(real64) :: u_cons(3), f_phys(3), rho1, u1, p1
    call prim_to_cons_flux(rho0, u0, p0, u_cons, f_phys)
    call cons_to_prim(u_cons, rho1, u1, p1)
    call assert_close(rho1, rho0, "roundtrip rho")
    call assert_close(u1,   u0,   "roundtrip u")
    call assert_close(p1,   p0,   "roundtrip p")
  end subroutine test_prim_cons_roundtrip

  !> flux_from_cons (flux from a conserved vector) must agree with the
  !> f_phys that prim_to_cons_flux produces from the matching primitives.
  !> This is the equivalence the Lax-Wendroff half-step relies on.
  subroutine test_flux_consistency()
    real(real64), parameter :: rho0 = 0.8_real64, u0 = -1.5_real64, p0 = 2.3_real64
    real(real64) :: u_cons(3), f_phys(3), f_from_cons(3)
    call prim_to_cons_flux(rho0, u0, p0, u_cons, f_phys)
    f_from_cons = flux_from_cons(u_cons)
    call assert_close(f_from_cons(1), f_phys(1), "flux consistency F(1)")
    call assert_close(f_from_cons(2), f_phys(2), "flux consistency F(2)")
    call assert_close(f_from_cons(3), f_phys(3), "flux consistency F(3)")
  end subroutine test_flux_consistency

end program check
