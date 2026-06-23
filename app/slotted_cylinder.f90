program run_slotted_cylinder
  use fluid_2d_models, only: slotted_cylinder
  implicit none
  integer :: val

  val = slotted_cylinder()
  if (val /= 1) then
    write(*, '(5X,A,/)') "Boo!!! Slotted cylinder 2d unsuccessful..."
    error stop 1
  end if
end program run_slotted_cylinder
