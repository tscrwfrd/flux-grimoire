module fluid_1d_models
  use iso_fortran_env, only: int32, real64, error_unit
  use fluid_forge
  implicit none

  private
  public :: square_wave, dam_break

contains

  !> Square wave function demonstrates how nicely the fct method transports
  !> and handles steep discontinuities pretty well, which can be useful in 
  !> problems that have shock waves.
  function square_wave() result(rtnvalue)
    real(real64), parameter :: dt = 0.000125
    real(real64), parameter :: dx = 1.0
    real(real64), parameter :: minnum = 1e-4
    real(real64), parameter :: dtdx = dt / dx
    integer(int32), parameter :: grid_size = 500
    real(real64), dimension(grid_size) :: qnty, velx, momx, temp_qnty, temp_momx
    integer(int32) :: rtnvalue, rc, funit, i, t

    open(action="write", file="./data/square_wave.csv", iostat=rc, newunit=funit,      &
      status="replace")

    if (rc /= 0) then
      write(error_unit, *) "Error opening a new file"
      rtnvalue = 0
      return
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
  end function square_wave

  !> This dam break simulation leverages the fct method and shows a weaknes in 
  !> in fct where the CFL condition must be set low to get good solutions. 
  !> Try setting it to different values...but beware of instabilities.
  function dam_break() result(rtnvalue)
    integer(int32) :: rtnvalue, t, i, funit, rc
    integer(int32), parameter :: nx = 504
    real(real64), parameter :: g = 9.8
    real(real64):: dt, dx, dtdx, cfl, max_wave_spd
    real(real64), dimension(nx) :: hgt, velx, momx, temp_hgt, &
      temp_mom

    momx(:) = 0.0
    velx(:) = 0.0
    ! initialize water levels
    hgt(1:250) = 20.0
    hgt(251:) = 0.5

    dx = 0.1 
    cfl = 0.15

    open(action="write", file="./data/dam_break.csv", iostat=rc, newunit=funit, &
      status="replace")

    if (rc /= 0) then
      write(error_unit, *) "Error opening a new file"
      rtnvalue = 0
      return
    end if

    do t = 1,1500
      ! max mave speed
      max_wave_spd = maxval(abs(velx) + sqrt(g*hgt))
      dt = cfl*dx/max_wave_spd
      dtdx = dt/dx
      
      ! advect water hgt
      ! ∂h/∂t + ∂(hu)/∂x
      call fct(hgt, velx, nx, dtdx, temp_hgt)

      ! advect  momentum
      ! ∂(hu)/∂t + ∂(hu²)/∂x
      momx = hgt * velx
      call fct(momx, velx, nx, dtdx, temp_mom)
      hgt = temp_hgt

      where(hgt > 0.001)
        velx = temp_mom / hgt
      elsewhere
        velx = 0.0
        hgt = 0.0
      end where
      ! leveraging operator splitting ∂(hu)/∂t + ∂(½gh²)/∂x
      ! central difference with updated velocity values from fct   
      ! remember fct has two padded cells on the boundaries
      where (hgt(3:nx-2) > 0.001)
        velx(3:nx-2) = velx(3:nx-2) - &
          (0.5*dtdx/hgt(3:nx-2))*(g*hgt(4:nx-1)**2 - &
          g*hgt(2:nx-3)**2)
      end where

      velx(1) = velx(3)
      velx(2) = velx(3)

      
      if (mod(t, 10) == 0) then 
        write(funit, '(*(g0.6,:,","))') hgt
        write(funit, '(*(g0.6,:,","))') velx
      end if
            
    end do 

    close(funit)
    rtnvalue = 1  

  end function dam_break

end module fluid_1d_models

