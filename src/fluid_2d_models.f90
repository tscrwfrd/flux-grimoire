module fluid_2d_models
  use iso_fortran_env, only: int32, real64, error_unit
  use fluid_forge, only: fct_2d
  implicit none

  private
  public :: slotted_cylinder

contains

  !> Zalesak's slotted-cylinder solid-body rotation test for the 2D FCT
  !> algorithm. A slotted disk on a unit-spacing 100x100 physical grid is
  !> advected through a rigid rotation velocity field. After one full
  !> revolution the disk should be back where it started; how cleanly the
  !> slot survives is the visual measure of the limiter's quality.
  !>
  !> Configuration follows the canonical setup:
  !>   - 100x100 physical cells + 2 ghost cells per side  -> nx = ny = 104
  !>   - cylinder center at (50, 75), radius 15
  !>   - slot of half-width 2.5 from y = 60 to y = 85
  !>   - rigid rotation about (50, 50) with omega = 1
  !>   - dt chosen so CFL stays well below the LW limit everywhere
  !>
  !> CSV layout: ny rows per snapshot (each row is one y-line of rho across
  !> all nx cells), starting with the initial condition.
  function slotted_cylinder() result(rtnvalue)
    integer(int32) :: rtnvalue, rc, funit, t, i, j
    integer(int32), parameter :: nx = 104, ny = 104
    integer(int32), parameter :: nsteps = 1250       ! ~2*pi / dt (one rotation)
    integer(int32), parameter :: snap_every = 50
    real(real64),   parameter :: dx = 1.0_real64, dy = 1.0_real64
    real(real64),   parameter :: dt = 0.005_real64
    real(real64),   parameter :: omega = 1.0_real64
    real(real64),   parameter :: xc = 50.0_real64, yc = 50.0_real64
    real(real64),   parameter :: cx = 50.0_real64, cy = 75.0_real64
    real(real64),   parameter :: r_cyl = 15.0_real64
    real(real64),   parameter :: slot_half = 2.5_real64
    real(real64),   parameter :: slot_top  = 85.0_real64
    real(real64),   parameter :: slot_bot  = 60.0_real64
    real(real64) :: x, y, dtdx, dtdy
    real(real64), dimension(nx, ny) :: rho, u, v, rho_new

    dtdx = dt / dx
    dtdy = dt / dy

    ! Cell-centered velocity field and slotted-disk initial condition.
    ! Physical coordinates are placed so the two ghost cells on each side
    ! sit at negative x/y and just past the 100-unit domain.
    rho = 0.0_real64
    do j = 1, ny
      do i = 1, nx
        x = (real(i, real64) - 2.5_real64) * dx
        y = (real(j, real64) - 2.5_real64) * dy
        u(i, j) = -omega * (y - yc)
        v(i, j) =  omega * (x - xc)
        if (sqrt((x - cx)**2 + (y - cy)**2) <= r_cyl) then
          if (.not. (abs(x - cx) <= slot_half .and. &
                     y >= slot_bot .and. y <= slot_top)) then
            rho(i, j) = 1.0_real64
          end if
        end if
      end do
    end do

    open(action="write", file="./data/slotted_cylinder.csv", iostat=rc, &
         newunit=funit, status="replace")
    if (rc /= 0) then
      write(error_unit, *) "Error opening a new file"
      rtnvalue = 0
      return
    end if

    do j = 1, ny
      write(funit, '(*(g0.6,:,","))') rho(:, j)
    end do

    rho_new = rho
    do t = 1, nsteps
      call fct_2d(rho, u, v, nx, ny, dtdx, dtdy, rho_new)
      rho = rho_new
      if (mod(t, snap_every) == 0) then
        do j = 1, ny
          write(funit, '(*(g0.6,:,","))') rho(:, j)
        end do
      end if
    end do

    close(funit)
    rtnvalue = 1
  end function slotted_cylinder

end module fluid_2d_models
