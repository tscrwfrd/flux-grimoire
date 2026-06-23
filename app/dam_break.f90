program run_dam_break
  use fluid_1d_models, only: dam_break
  implicit none
  integer :: val

  val = dam_break()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Dam break 1d unsuccessful..."
    error stop 1
  end if
end program run_dam_break
