program main
  use fluid_forge, only: say_hello, fct
  use run1d, only: sim1d
  implicit none
  integer :: val

  call say_hello()
  val = sim1d()

end program main
