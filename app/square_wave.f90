program run_square_wave
  use fluid_1d_models, only: square_wave
  implicit none
  integer :: val

  val = square_wave()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! square wave 1d unsuccessful..."
    error stop 1
  end if
end program run_square_wave
