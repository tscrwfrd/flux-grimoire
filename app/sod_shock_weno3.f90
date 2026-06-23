program run_sod_shock_weno3
  use fluid_1d_models, only: sod_shock_weno3
  implicit none
  integer :: val

  val = sod_shock_weno3()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock WENO3 1d unsuccessful..."
    error stop 1
  end if
end program run_sod_shock_weno3
