module run1d
  use iso_fortran_env, only: int32, real64, error_unit
  use fluid_forge
  implicit none

  private
  public :: sim1d

contains

  function sim1d() result(rtnvalue)
    real(real64), parameter :: dt = 0.000125
    real(real64), parameter :: dx = 1.0
    real(real64), parameter :: minnum = 1e-4
    real(real64), parameter :: dtdx = dt / dx
    integer(int32), parameter :: grid_size = 500
    real(real64), dimension(grid_size) :: qnty, velx, momx, temp_qnty, temp_momx
    integer(int32) :: rtnvalue, rc, funit, i, t

    open(action="write", file="fct_out.csv", iostat=rc, newunit=funit,      &
      status="replace")

    if (rc /= 0) then
      write(error_unit, *) "Error opening a new file"
      stop
    end if 


    ! initialize square wave
    qnty(50:100) = 250.0
    velx(50:100) = 1000.0

    do t = 1, 5000
      ! transport 
      call fct(qnty, velx, grid_size, dtdx, temp_qnty)
      momx = qnty * velx

      call fct(momx, velx, grid_size, dtdx,temp_momx)
      qnty = temp_qnty
      momx = temp_momx

      where (qnty > minnum)
        velx = momx / qnty
      elsewhere
        velx = 0.0
        qnty = 0.0
      end where

      if (mod(t, 10) == 0) then 
        write(funit, '(*(g0.6,:,","))') qnty
        write(funit, '(*(g0.6,:,","))') velx
      end if
    
    end do 
    close(funit)
    rtnvalue = 1
  end function sim1d

end module run1d

