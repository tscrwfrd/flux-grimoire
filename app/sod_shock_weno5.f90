program run_sod_shock_weno5
  use fluid_1d_models, only: sod_shock_weno5
  implicit none
  integer :: val

  val = sod_shock_weno5()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock WENO5 1d unsuccessful..."
    error stop 1
  end if
end program run_sod_shock_weno5
