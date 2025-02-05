module run1d_dam_break
  use iso_fortran_env, only: int32, real64, error_unit
  use fluid_forge
  implicit none

  private
  public:: dam1d

contains 

  function dam1d() result(rtnvalue)
    integer(int32) :: rtnvalue, t, i, funit, rc
    integer(int32), parameter :: grid_size = 504
    real(real64), parameter :: g = 9.8
    real(real64):: dt, dx, dtdx, cfl, max_wave_spd
    real(real64), dimension(grid_size) :: height, velx, momx, temp_height, &
      temp_mom

    ! initialize
    momx(:) = 0.0
    velx(:) = 0.0
    height(1:250) = 10.0
    height(251:) = 0.5

    open(action="write", file="dam_out.csv", iostat=rc, newunit=funit, &
      status="replace")

    if (rc /= 0) then
      write(error_unit, *) "Error opening a new file"
      stop
    end if

    dt = 0.0001
    dx = 0.1 
    dtdx = dt / dx
    cfl = 0.2

    do t = 1,1500
      max_wave_spd = maxval(abs(velx) + sqrt(g*height))

      ! if (max_vel > 0.01) then
      dt = cfl*dx/max_wave_spd
      ! print*, cdt
      ! dt = 0.0001
      ! else
        ! dt = 0.0001
      ! end if

      dtdx = dt/dx
      ! advect water height
      call fct(height, velx, grid_size, dtdx, temp_height)
      momx = height * velx

      ! advect  momentum
      call fct(momx, velx, grid_size, dtdx, temp_mom)
      height = temp_height
      momx = temp_mom

      where(height > 0.001)
        velx = momx / height
      elsewhere
        velx = 0.0
        height = 0.0
      end where

      ! dtdx = dt/dx
      ! apply pressure gradient
      ! where(height > 0.001)
      
      where (height(3:502) > 0.001)
        temp_mom(3:502) = velx(3:502) - &
          (0.5*dtdx/height(3:502))*(g*height(4:503)**2 - g*height(2:501)**2)
      ! elsewhere
      !   temp_mom = 0.0
      end where
      ! temp_mom(3:502) = velx(3:502) - &
      !   (0.5*dtdx/height(3:502))*(g*height(4:503)**2 - g*height(2:501)**2)

      velx = temp_mom
      

      if (mod(t, 10) == 0) then 
        write(funit, '(*(g0.6,:,","))') height
        write(funit, '(*(g0.6,:,","))') velx
      end if
      
            
    end do 

    close(funit)
    rtnvalue = 1  
    


  end function dam1d

  
end module run1d_dam_break
