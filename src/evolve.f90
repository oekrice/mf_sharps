!*******************************************************************************

MODULE evolve
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    USE boundary
    USE pressure
    USE netcdf

    !USE output
    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE timestep()

    !New timestep with a quasi-implicit bit. Using a predictor for the magnfield of the new bit but not the old
    !if (n < 17) then

    CALL calculate_magnetic()

    !CALL check_solenoidal()

    CALL calculate_jp()

    CALL calculate_current()

    !print*, 'j', maxval(abs(jx)),  maxval(abs(jy)),  maxval(abs(jz))
    CALL j_to_gridpts()
    CALL b_to_gridpts()

    CALL calculate_velocity()

    !CALL calculate_pressure()
    CALL calculate_electric()
    !end if

    CALL MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.

END SUBROUTINE timestep


SUBROUTINE calculate_magnetic()

    IMPLICIT NONE

    !Determine the interior points from the vector potential.
    !No boundary conditions here (do that in a bit)

    ! INTERIOR POINTS (DON'T USE INFORMATION FROM A THAT DOESN'T EXIST)

    bx(0:nx, 1:ny,1:nz) = (az(0:nx,1:ny,1:nz) - az(0:nx, 0:ny-1,1:nz))/dy - (ay(0:nx,1:ny,1:nz) - ay(0:nx,1:ny,0:nz-1))/dz

    by(1:nx, 0:ny,1:nz) = (ax(1:nx,0:ny,1:nz) - ax(1:nx, 0:ny,0:nz-1))/dz - (az(1:nx,0:ny,1:nz) - az(0:nx-1,0:ny,1:nz))/dx

    bz(1:nx, 1:ny,0:nz) = (ay(1:nx,1:ny,0:nz) - ay(0:nx-1, 1:ny,0:nz))/dx - (ax(1:nx,1:ny,0:nz) - ax(1:nx,0:ny-1,0:nz))/dy

    CALL bfield_mpi
    CALL magnetic_boundary

    if (n == 0) bz_surf_reference(0:nx+1,0:ny+1) = bz(0:nx+1,0:ny+1,0)  !Save reference lower boundary field to stop annoying instabilities due to lack of upwinding

END SUBROUTINE calculate_magnetic

SUBROUTINE calculate_current()

    IMPLICIT NONE
    !Determine the current from the magnetic field (after boundary conditions etc.)

    jx(0:nx+1, 0:ny,0:nz) = (bz(0:nx+1,1:ny+1,0:nz) - bz(0:nx+1, 0:ny,0:nz))/dy - (by(0:nx+1,0:ny,1:nz+1) - by(0:nx+1,0:ny,0:nz))/dz

    jy(0:nx, 0:ny+1,0:nz) = (bx(0:nx,0:ny+1,1:nz+1) - bx(0:nx, 0:ny+1,0:nz))/dz - (bz(1:nx+1,0:ny+1,0:nz) - bz(0:nx,0:ny+1,0:nz))/dx

    jz(0:nx, 0:ny,0:nz+1) = (by(1:nx+1,0:ny,0:nz+1) - by(0:nx, 0:ny,0:nz+1))/dx - (bx(0:nx,1:ny+1,0:nz+1) - bx(0:nx,0:ny,0:nz+1))/dy

    jx(0:nx+1, 0:ny,0:nz) = jx(0:nx+1, 0:ny,0:nz) - jpx(0:nx+1, 0:ny,0:nz)
    jy(0:nx, 0:ny+1,0:nz) = jy(0:nx, 0:ny+1,0:nz) - jpy(0:nx, 0:ny+1,0:nz)

END SUBROUTINE calculate_current

SUBROUTINE j_to_gridpts
    !Averages the current field to raw grid points
    !Should only need to average in one direction
    IMPLICIT NONE
    jx1(0:nx,0:ny,0:nz) = 0.5_num*(jx(1:nx+1,0:ny,0:nz) + jx(0:nx,0:ny,0:nz))
    jy1(0:nx,0:ny,0:nz) = 0.5_num*(jy(0:nx,1:ny+1,0:nz) + jy(0:nx,0:ny,0:nz))
    jz1(0:nx,0:ny,0:nz) = 0.5_num*(jz(0:nx,0:ny,1:nz+1) + jz(0:nx,0:ny,0:nz))

END SUBROUTINE j_to_gridpts

SUBROUTINE b_to_gridpts
    !Averages the magnetic field to raw grid points
    !Need to average in two dimensions
    IMPLICIT NONE
    bx1(0:nx,0:ny,0:nz) = 0.25_num*(bx(0:nx,0:ny,0:nz) + bx(0:nx,1:ny+1,0:nz) + bx(0:nx,0:ny,1:nz+1) + bx(0:nx,1:ny+1,1:nz+1))
    by1(0:nx,0:ny,0:nz) = 0.25_num*(by(0:nx,0:ny,0:nz) + by(1:nx+1,0:ny,0:nz) + by(0:nx,0:ny,1:nz+1) + by(1:nx+1,0:ny,1:nz+1))
    bz1(0:nx,0:ny,0:nz) = 0.25_num*(bz(0:nx,0:ny,0:nz) + bz(1:nx+1,0:ny,0:nz) + bz(0:nx,1:ny+1,0:nz) + bz(1:nx+1,1:ny+1,0:nz))

END SUBROUTINE b_to_gridpts

SUBROUTINE calculate_velocity
    !Calculates the magnetofrictional velocity
    IMPLICIT NONE
    real(num), dimension(:,:,:):: b2_lim(0:nx,0:ny,0:nz)
    real(num), dimension(:,:,:):: bx_lim(0:nx,0:ny,0:nz), by_lim(0:nx,0:ny,0:nz), bz_lim(0:nx,0:ny,0:nz)
    real(num), dimension(:,:,:):: jx_lim(0:nx,0:ny,0:nz), jy_lim(0:nx,0:ny,0:nz), jz_lim(0:nx,0:ny,0:nz)
    real(num):: min_value

    !Stop underflow errors here (I know, it's a bit of a shame)
    min_value = 1D-8
    bx_lim = merge(0.0_num, bx1, abs(bx1) < min_value)
    by_lim = merge(0.0_num, by1, abs(by1) < min_value)
    bz_lim = merge(0.0_num, bz1, abs(bz1) < min_value)

    jx_lim = merge(0.0_num, jx1, abs(jx1) < min_value)
    jy_lim = merge(0.0_num, jy1, abs(jy1) < min_value)
    jz_lim = merge(0.0_num, jz1, abs(jz1) < min_value)

    b2 = bx_lim**2 + by_lim**2 + bz_lim**2 !B squared

    nu(:,:,:) = nu0

    if (abs(mf_delta) < 1e-10) then !No softening
        soft = b2
    else !Softening.
        !Get rid of underflow error in the exponential using b2_lim, hopefully
        b2_lim = merge(10_num*mf_delta, b2, b2 > 10_num*mf_delta)
        soft = b2 + mf_delta*exp(-b2_lim/mf_delta)
    end if

    vx = nu*(jy_lim*bz_lim - jz_lim*by_lim)/soft
    vy = nu*(jz_lim*bx_lim - jx_lim*bz_lim)/soft
    vz = nu*(jx_lim*by_lim - jy_lim*bx_lim)/soft

    vx = merge(vx, max_velocity*sign(1.0_num, vx), abs(vx) < max_velocity)
    vy = merge(vy, max_velocity*sign(1.0_num, vy), abs(vy) < max_velocity)
    vz = merge(vz, max_velocity*sign(1.0_num, vz), abs(vz) < max_velocity)

    !print*, max_velocity, minval(abs(soft)), maxval(abs(vx)), maxval(abs(vy)), maxval(abs(vz))

    if (z_down < 0) then
        vx(:,:,0) = 0.0_num; vy(:,:,0) = 0.0_num; vz(:,:,0) = 0.0_num
    end if

END SUBROUTINE calculate_velocity

SUBROUTINE add_boundary_flows()

    !Adds twisting directly onto the electric field. Only affects ez (for the jet model, at least)
    IMPLICIT NONE

    if (z_down < 0 .and. shearfact > 1d-6) then

      vx(0:nx,0:ny,0) = shearfact*surf_vx(0:nx,0:ny)
      vy(0:nx,0:ny,0) = shearfact*surf_vy(0:nx,0:ny)
      vz(0:nx,0:ny,0) = shearfact*surf_vz(0:nx,0:ny)

    end if

END SUBROUTINE add_boundary_flows


SUBROUTINE check_solenoidal()
    !Checks the solenoidal condition by calculating the divergence of all raw grid cells
    IMPLICIT NONE

    real(num), dimension(:,:,:):: div(0:nx+1,0:ny+1,0:nz+1)
    div = 0.0_num
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dx*dy*(bz(0:nx+1,0:ny+1,-1:nz) - bz(0:nx+1,0:ny+1,0:nz+1))
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dy*dz*(bx(-1:nx,0:ny+1,0:nz+1) - bx(0:nx+1,0:ny+1,0:nz+1))
    div(0:nx+1,0:ny+1,0:nz+1) = div(0:nx+1,0:ny+1,0:nz+1) + dx*dz*(by(0:nx+1,-1:ny,0:nz+1) - by(0:nx+1,0:ny+1,0:nz+1))

    print*, 'Max divergence', maxval(abs(div(1:nx,1:ny,1:nz)))

END SUBROUTINE check_solenoidal

SUBROUTINE calculate_electric()

    !Calculates the electric field - resistivity, magnetofriction and boundary effects
    IMPLICIT NONE
    real(num), dimension(:,:,:):: jx_lim(0:nx,0:ny,0:nz), jy_lim(0:nx,0:ny,0:nz), jz_lim(0:nx,0:ny,0:nz)
    real(num):: min_value

    min_value = 1D-8

    ex = 0.0; ey = 0.0; ez = 0.0

    jx_lim = jx!merge(0.0_num, jx, abs(jx) < min_value)
    jy_lim = jy!merge(0.0_num, jy, abs(jy) < min_value)
    jz_lim = jz!merge(0.0_num, jz, abs(jz) < min_value)

    if (eta > 0) then
        !Determine the current from the magnetic field (after boundary conditions etc.)
        ex(1:nx, 0:ny,0:nz) = ex(1:nx, 0:ny,0:nz) + eta*(jx_lim(1:nx, 0:ny,0:nz))! - jpx(1:nx, 0:ny,0:nz))
        ey(0:nx, 1:ny,0:nz) = ey(0:nx, 1:ny,0:nz) + eta*(jy_lim(0:nx, 1:ny,0:nz))! - jpy(0:nx, 1:ny,0:nz))
        ez(0:nx, 0:ny,1:nz) = ez(0:nx, 0:ny,1:nz) + eta*jz_lim(0:nx, 0:ny,1:nz)
    end if

    !Add shearing (if necessary) directly onto this (averaged) field
    CALL add_boundary_flows()

    ex1 = vz*by1 - vy*bz1
    ey1 = vx*bz1 - vz*bx1
    ez1 = vy*bx1 - vx*by1

    !Average to Ribs (interior only):
    ex(1:nx,0:ny,0:nz) = ex(1:nx,0:ny,0:nz)  + 0.5_num*(ex1(0:nx-1,0:ny,0:nz) + ex1(1:nx,0:ny,0:nz))
    ey(0:nx,1:ny,0:nz) = ey(0:nx,1:ny,0:nz)  + 0.5_num*(ey1(0:nx,0:ny-1,0:nz) + ey1(0:nx,1:ny,0:nz))
    ez(0:nx,0:ny,1:nz) = ez(0:nx,0:ny,1:nz)  + 0.5_num*(ez1(0:nx,0:ny,0:nz-1) + ez1(0:nx,0:ny,1:nz))

    !Add outflow (if necessary) directly onto this field
    if (voutfact > 0) then
    ex(1:nx,0:ny,0:nz) = ex(1:nx,0:ny,0:nz) + voutx(1:nx,0:ny,0:nz)*by(1:nx,0:ny,0:nz)
    ey(0:nx,1:ny,0:nz) = ey(0:nx,1:ny,0:nz) - vouty(0:nx,1:ny,0:nz)*bx(0:nx,1:ny,0:nz)
    end if
    
    !Add electric field loaded in from elsewhere
    if (z_rank == 0) then
        ex(1:nx,0:ny,0) = surf_ex(1:nx,0:ny)
        ey(0:nx,1:ny,0) = surf_ey(0:nx,1:ny)
    end if


END SUBROUTINE calculate_electric

SUBROUTINE import_surface_electric(flow_number, dt_fact)

    !Imports the velocity field from the imported DAVE magnetogram
    IMPLICIT NONE

    INTEGER:: flow_number

    CHARACTER(LEN =64):: electric_filename
    CHARACTER(LEN = 4):: flow_id
    CHARACTER(LEN = 4):: run_id

    INTEGER:: ncid, vid
    REAL(num):: dt_fact

    if (flow_number < 499) then
        write (flow_id,'(I4.4)') flow_number
        write (run_id,'(I3.3)') init_number

        electric_filename = trim("./efields/"//trim(run_id)//'/'//trim(flow_id)//'.nc')

        call try(nf90_open(trim(electric_filename), nf90_nowrite, ncid))

        call try(nf90_inq_varid(ncid, 'ex', vid))
        call try(nf90_get_var(ncid, vid, surf_ex(0:nx+1,0:ny), &
        start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+2,ny+1/)))

        call try(nf90_inq_varid(ncid, 'ey', vid))
        call try(nf90_get_var(ncid, vid, surf_ey(0:nx,0:ny+1), &
        start = (/x_rank*nx+1,y_rank*ny+1/),count = (/nx+1,ny+2/)))

        call try(nf90_close(ncid))

        surf_ex = surf_ex*dt_fact
        surf_ey = surf_ey*dt_fact

    else
        surf_ex = 0.0_num
        surf_ey = 0.0_num
    end if

END SUBROUTINE import_surface_electric

subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        call mpi_abort(comm, ierr)
    end if

end subroutine try



!*******************************************************************************
END MODULE evolve
!*******************************************************************************
