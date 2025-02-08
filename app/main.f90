program main
  use fluid_forge, only: say_hello
  use fluid_1d_models, only: dam_break, square_wave, sod_shock
  implicit none
  integer :: val

  call say_hello()

  val = square_wave()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! square wave 1d unsuccessful..."
  end if
  
  val = dam_break()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Dam break 1d unsuccessful..."
  end if

  val = sod_shock()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock 1d unsuccessful..."
  end if

end program main
