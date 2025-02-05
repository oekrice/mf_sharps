!*******************************************************************************

MODULE pressure
!*******************************************************************************
! Functions for the new pressure part
! Gives a function to add on to the velocity term, which is averaged to grid POINTS

!*******************************************************************************
    USE shared_data

    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE pressure_function()

    !Input the pressure function f(y) here - note y is the vertical coordinate like LARE
    !Should this be positive or negative?!?!
    implicit none
    integer:: k

    do k = 0, nz+1
        if (zstar < 0.01) then
            if (k == 0 .and. proc_num == -1 .and. mag_min == 0) print*, 'Base case, no pressure function'
            fz(:,:,k) = 0.0_num
        else
            if (k == 0 .and. proc_num == 0) then
                print*, '___________________'
                print*, 'Pressure parameters'
                print*, 'a = ', a
                print*, 'b = ', b
                print*, 'zstar = ', zstar
                print*, 'deltaz = ', deltaz
            end if
            if (decay_type == 0) then
                if (k == 0 .and. proc_num == 0) print*, 'Base case, no pressure function'
                    fz(:,:,k) = 0.0_num
            else if (decay_type == 1) then  !Exponential
                if (k == 0 .and. proc_num == 0) print*, 'Exponential Pressure Function'
                    fz(:,:,k) = a*exp(-zs(k)/b)
            else if (decay_type == 2) then !Smooth tanh
                if (k == 0 .and. proc_num == 0) print*, 'Smooth Tanh Pressure Function'
                    fz(:,:,k) = a*(1.0_num - b*tanh((zs(k)-zstar)/deltaz))
            else !Sharp tanh
                if (k == 0 .and. proc_num == 0) print*, 'Sharp Tanh Pressure Function'
                    fz(:,:,k) = a*(1.0_num - b*tanh((zs(k)-zstar)/deltaz))
            end if
        end if
    end do


END SUBROUTINE pressure_function

SUBROUTINE calculate_jp()
    !Calculates jp - the 'pressure current'
    implicit none

    jpx(0:nx+1, 0:ny,0:nz) = (fz(0:nx+1,1:ny+1,0:nz)*bz(0:nx+1,1:ny+1,0:nz) - fz(0:nx+1,0:ny,0:nz)*bz(0:nx+1, 0:ny,0:nz))/dy
    jpy(0:nx, 0:ny+1,0:nz) = -(fz(1:nx+1,0:ny+1,0:nz)*bz(1:nx+1,0:ny+1,0:nz) - fz(0:nx,0:ny+1,0:nz)*bz(0:nx,0:ny+1,0:nz))/dx

    !Stop instabilities on the surface by setting the 'pressure diffusion' to zero here
    if (z_down < 0) then
        jpx(:,:,0) = 0.0_num; jpy(:,:,0) = 0.0_num
    end if

    jpx1(0:nx,0:ny,0:nz) = 0.5_num*(jpx(1:nx+1,0:ny,0:nz) + jpx(0:nx,0:ny,0:nz))
    jpy1(0:nx,0:ny,0:nz) = 0.5_num*(jpy(0:nx,1:ny+1,0:nz) + jpy(0:nx,0:ny,0:nz))

END SUBROUTINE calculate_jp

SUBROUTINE calculate_pressure()
    !Does the cross product with b, averages to gridpoints and does the softening as for the velocity
    !It appears that it might be necessary to do this step rather than just adding on to the current.
    !There must have been a reason why I did it in the 2D cases...

    !Can get nu directly from the velocity calculation, so don't need to do this twice. But this function must come after calculate_velocity.

    implicit none

    REAL(num), DIMENSION(0:nx,0:ny,0:nz):: vpx, vpy, vpz

    !Extra velocity due to the `pressure current'.
    vpx = nu*(jpy1*bz1 - 0.0_num )/soft
    vpy = nu*(0.0_num  - jpx1*bz1)/soft
    vpz = nu*(jpx1*by1 - jpy1*bx1)/soft

    !Subtract from the already-calculated velocity
    vx(0:nx,0:ny,0:nz) = vx(0:nx,0:ny,0:nz) - vpx(0:nx, 0:ny,0:nz)
    vy(0:nx,0:ny,0:nz) = vy(0:nx,0:ny,0:nz) - vpy(0:nx, 0:ny,0:nz)
    vz(0:nx,0:ny,0:nz) = vz(0:nx,0:ny,0:nz) - vpz(0:nx, 0:ny,0:nz)

    if (z_down < 0) then
        vpx(:,:,0) = 0.0_num; vpy(:,:,0) = 0.0_num; vpz(:,:,0) = 0.0_num
    end if

END SUBROUTINE calculate_pressure





!*******************************************************************************
END MODULE pressure
!*******************************************************************************
