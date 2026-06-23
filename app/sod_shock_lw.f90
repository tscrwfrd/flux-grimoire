program run_sod_shock_lw
  use fluid_1d_models, only: sod_shock_lw
  implicit none
  integer :: val

  val = sod_shock_lw()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock LW 1d unsuccessful..."
    error stop 1
  end if
end program run_sod_shock_lw
