program main
  use fluid_forge, only: say_hello
  use fluid_1d_models, only: dam_break, square_wave, sod_shock_lw, sod_shock_roe, &
    sod_shock_weno3, sod_shock_weno5
  use fluid_2d_models, only: slotted_cylinder
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

  val = sod_shock_lw()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock LW 1d unsuccessful..."
  end if

  val = sod_shock_roe()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock Roe 1d unsuccessful..."
  end if

  val = sod_shock_weno3()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock WENO3 1d unsuccessful..."
  end if

  val = sod_shock_weno5()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Sod shock WENO5 1d unsuccessful..."
  end if

  val = slotted_cylinder()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Slotted cylinder 2d unsuccessful..."
  end if

end program main
