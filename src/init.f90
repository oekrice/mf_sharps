!*******************************************************************************
MODULE init
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    USE netcdf
    USE pressure

    !USE output
    IMPLICIT NONE

!*******************************************************************************

CONTAINS

SUBROUTINE initialise()

    CALL read_parameters()

    CALL start_mpi()

    CALL MPI_BARRIER(comm, ierr)
    if (proc_num == -1) then
        print*, 'Run Number', int(run_number)
        print*, 'Outflow', voutfact
        print*, 'Shearing', shearfact
        print*, 'eta', eta
        print*, 'nu0', nu0
        print*, 'eta0', eta0
        print*, 'Data directory', data_directory
    end if

    CALL MPI_create_types

    CALL allocate_arrays()

    CALL establish_grid()

    CALL calculate_timestep()

    CALL set_outflow()

    CALL set_shearing()

    CALL pressure_function()

END SUBROUTINE initialise

SUBROUTINE read_parameters()
    !Reads in the parameters from the text file variables.txt
    REAL(num), DIMENSION(40):: variables(0:39)
    CHARACTER(LEN=64):: input_value
    CHARACTER(LEN=4):: run_id

    CHARACTER(LEN=64):: parameter_filename, magtimes_filename

    call get_command_argument(1, input_value)
    read(unit=input_value,fmt=*) run_number

    !Make magnetogram filename
    write (run_id,'(I3.3)') int(run_number)

    parameter_filename = trim('./parameters/variables'//trim(run_id)//'.txt')
    magtimes_filename = trim('./parameters/magtimes'//trim(run_id)//'.txt')

    OPEN(1, FILE = parameter_filename)
    READ(1, *) variables
    CLOSE(1)

    run_number = variables(0)
    nx_global = int(variables(1))
    ny_global = int(variables(2))
    nz_global = int(variables(3))
    tmax = variables(4)

    nplots = variables(5)
    ndiags = int(variables(6))

    voutfact = variables(7)
    shearfact = variables(8)
    eta = variables(9)
    nu0 = variables(10)
    eta0 = variables(11)

    x0_global = variables(12)
    x1_global = variables(13)
    y0_global = variables(14)
    y1_global = variables(15)
    z0_global = variables(16)
    z1_global = variables(17)

    a = variables(18)
    b = variables(19)
    deltaz = variables(20)
    zstar = variables(21)

    hamilton_flag = int(variables(22))
    decay_type = int(variables(23))

    continue_time = variables(24)

    nmags = int(variables(25))

    init_number = int(variables(27))

    mag_min = int(variables(29))
    mag_max = int(variables(30))

    !Read in magnetogram times (need to be able to vary the times here)
    allocate(mag_times(0:nmags-1))

    OPEN(1, FILE = magtimes_filename)
    READ(1, *) mag_times
    CLOSE(1)

END SUBROUTINE read_parameters

SUBROUTINE allocate_arrays()
    ! - Estalish the grid (not particularly difficult for a cartesian box). Try to use the same numbering system as LARE.
    ! - index 0 is always the cell on the boundary, or just OUTSIDE it. c cells are numbered from -1, s cells from -2
    ! - Arrays are local unless otherwise declared. Not sure if any global arrays are required, but perhaps

    allocate(xs(-1:nx+1)); allocate(ys(-1:ny+1)); allocate(zs(-1:nz+1))
    allocate(xc(0:nx+1)); allocate(yc(0:ny+1)); allocate(zc(0:nz+1))

    allocate(xs_global(-1:nx_global+1)); allocate(ys_global(-1:ny_global+1)); allocate(zs_global(-1:nz_global+1))
    allocate(xc_global(0:nx_global+1)); allocate(yc_global(0:ny_global+1)); allocate(zc_global(0:nz_global+1))


    !Allocate global arrays (things that could conceivably be plotted afterwards, I suppose)
    allocate(ax(0:nx+1,-1:ny+1,-1:nz+1)); allocate(ay(-1:nx+1,0:ny+1,-1:nz+1)); allocate(az(-1:nx+1,-1:ny+1,0:nz+1))
    allocate(bx(-1:nx+1,0:ny+1,0:nz+1)); allocate(by(0:nx+1,-1:ny+1,0:nz+1)); allocate(bz(0:nx+1,0:ny+1,-1:nz+1))
    allocate(jx(0:nx+1,-1:ny+1,-1:nz+1)); allocate(jy(-1:nx+1,0:ny+1,-1:nz+1)); allocate(jz(-1:nx+1,-1:ny+1,0:nz+1))
    allocate(ex(0:nx+1,-1:ny+1,-1:nz+1)); allocate(ey(-1:nx+1,0:ny+1,-1:nz+1)); allocate(ez(-1:nx+1,-1:ny+1,0:nz+1))

    !Fields averaged to gridpoints
    allocate(jx1(0:nx,0:ny,0:nz)); allocate(jy1(0:nx,0:ny,0:nz)); allocate(jz1(0:nx,0:ny,0:nz))
    allocate(bx1(0:nx,0:ny,0:nz)); allocate(by1(0:nx,0:ny,0:nz)); allocate(bz1(0:nx,0:ny,0:nz))
    allocate(ex1(0:nx,0:ny,0:nz)); allocate(ey1(0:nx,0:ny,0:nz)); allocate(ez1(0:nx,0:ny,0:nz))

    allocate(vx(0:nx,0:ny,0:nz)); allocate(vy(0:nx,0:ny,0:nz)); allocate(vz(0:nx,0:ny,0:nz))
    allocate(vx_surf(0:nx,0:ny)); allocate(vy_surf(0:nx,0:ny))

    allocate(lf(0:nx,0:ny,0:nz)); allocate(b2(0:nx,0:ny,0:nz)); allocate(nu(0:nx,0:ny,0:nz))
    allocate(soft(0:nx,0:ny,0:nz))

    allocate(bz_surf_reference(0:nx+1,0:ny+1))

    !Pressure arrays
    allocate(fz(0:nx+1,0:ny+1,-1:nz+1)) !Aligned with bz (I'm pretty sure)
    allocate(jpx(0:nx+1,-1:ny+1,-1:nz+1)); allocate(jpy(-1:nx+1,0:ny+1,-1:nz+1))
    allocate(jpx1(0:nx,0:ny,0:nz)); allocate(jpy1(0:nx,0:ny,0:nz))

    allocate(surf_vx(0:nx,0:ny)); allocate(surf_vy(0:nx,0:ny)); allocate(surf_vz(0:nx,0:ny))
    allocate(surf_ex(0:nx+1,-1:ny+1)); allocate(surf_ey(-1:nx+1,0:ny+1))

    allocate(surf_bx0(-1:nx+1,0:ny+1)); allocate(surf_bx1(-1:nx+1,0:ny+1))
    allocate(surf_by0(0:nx+1,-1:ny+1)); allocate(surf_by1(0:nx+1,-1:ny+1))
    allocate(surf_bx(-1:nx+1,0:ny+1)); allocate(surf_by(0:nx+1,-1:ny+1))

    bx = 0.0_num; by = 0.0_num; bz = 0.0_num

END SUBROUTINE allocate_arrays

SUBROUTINE establish_grid()

    ! Imports the initial condition from the inits file, IF mag_min is 0. Otherwise read in a previously-saved snapshot
    ! Will try to be smart with this, such that the vector potential is only read in in individual processes

    CHARACTER(LEN =64):: init_filename
    INTEGER:: ncid, vid
    CHARACTER(LEN = 4):: mag_id, run_id, init_id

    write (mag_id,'(I4.4)') mag_min
    write (run_id,'(I3.3)') int(run_number)
    write (init_id,'(I3.3)') int(init_number)

    if (mag_min == 0) then
        init_filename = trim('./inits/init'//trim(init_id)//'.nc')
    else
        init_filename = trim(trim(data_directory_root)//trim(run_id)//'/'//trim(mag_id)//'.nc')
    end if


    if (proc_num == -1) print*, run_id, mag_id, 'Init filename:  ', init_filename
    call try(nf90_open(trim(init_filename), nf90_nowrite, ncid))

    call try(nf90_inq_varid(ncid, 'ax', vid))
    call try(nf90_get_var(ncid, vid, ax(1:nx,0:ny,0:nz), &
    start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx,ny+1,nz+1/)))

    call try(nf90_inq_varid(ncid, 'ay', vid))
    call try(nf90_get_var(ncid, vid, ay(0:nx,1:ny,0:nz), &
    start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny,nz+1/)))

    call try(nf90_inq_varid(ncid, 'az', vid))
    call try(nf90_get_var(ncid, vid, az(0:nx,0:ny,1:nz), &
    start = (/x_rank*nx+1,y_rank*ny+1,z_rank*nz+1/),count = (/nx+1,ny+1,nz/)))


    call try(nf90_inq_varid(ncid, 'xs', vid))
    call try(nf90_get_var(ncid, vid, xs(0:nx), &
    start = (/x_rank*nx+1/),count = (/nx+1/)))

    call try(nf90_inq_varid(ncid, 'ys', vid))
    call try(nf90_get_var(ncid, vid, ys(0:ny), &
    start = (/y_rank*ny+1/),count = (/ny+1/)))

    call try(nf90_inq_varid(ncid, 'zs', vid))
    call try(nf90_get_var(ncid, vid, zs(0:nz), &
    start = (/z_rank*nz+1/),count = (/nz+1/)))


    xs(-1) = 2*xs(0) - xs(1); xs(nx+1) = 2*xs(nx) - xs(nx-1)
    ys(-1) = 2*ys(0) - ys(1); ys(ny+1) = 2*ys(ny) - ys(ny-1)
    zs(-1) = 2*zs(0) - zs(1); zs(nz+1) = 2*zs(nz) - zs(nz-1)

    dx = sum(xs(0:nx+1) - xs(-1:nx))/(nx+2)
    dy = sum(ys(0:ny+1) - ys(-1:ny))/(ny+2)
    dz = sum(zs(0:nz+1) - zs(-1:nz))/(nz+2)

    xc = 0.5_num*(xs(-1:nx) + xs(0:nx+1))
    yc = 0.5_num*(ys(-1:ny) + ys(0:ny+1))
    zc = 0.5_num*(zs(-1:nz) + zs(0:nz+1))

    x0 = xs(0); x1 = xs(nx)
    y0 = ys(0); y1 = ys(ny)
    z0 = zs(0); z1 = zs(nz)

    call try(nf90_inq_varid(ncid, 'xs', vid))
    call try(nf90_get_var(ncid, vid, xs_global(0:nx_global), &
    start = (/1/),count = (/nx_global+1/)))

    call try(nf90_inq_varid(ncid, 'ys', vid))
    call try(nf90_get_var(ncid, vid, ys_global(0:ny_global), &
    start = (/1/),count = (/ny_global+1/)))

    call try(nf90_inq_varid(ncid, 'zs', vid))
    call try(nf90_get_var(ncid, vid, zs_global(0:nz_global), &
    start = (/1/),count = (/nz_global+1/)))

    call try(nf90_close(ncid))

    xs_global(-1) = 2*xs_global(0) - xs_global(1)
    xs_global(nx_global+1) = 2*xs_global(nx_global) - xs_global(nx_global-1)
    ys_global(-1) = 2*ys_global(0) - ys_global(1)
    ys_global(ny_global+1) = 2*ys_global(ny_global) - ys_global(ny_global-1)
    zs_global(-1) = 2*zs_global(0) - zs_global(1)
    zs_global(nz_global+1) = 2*zs_global(nz_global) - zs_global(nz_global-1)

    xc_global = 0.5_num*(xs_global(-1:nx_global) + xs_global(0:nx_global+1))
    yc_global = 0.5_num*(ys_global(-1:ny_global) + ys_global(0:ny_global+1))
    zc_global = 0.5_num*(zs_global(-1:nz_global) + zs_global(0:nz_global+1))

    volume = (x1-x0)*(y1-y0)*(z1-z0)
    volume_global = (x1_global-x0_global)*(y1_global-y0_global)*(z1_global-z0_global)

END SUBROUTINE establish_grid

subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        call mpi_abort(comm, ierr)
    end if

end subroutine try

SUBROUTINE calculate_timestep()
    ! Calculates the timestep based on the parameters and cfl condition. Then adjusts so there are an integer number of steps in between snapshots, as that's neater
    REAL(num):: dt_ideal, nt_ft
    REAL(num):: plot_dt

    dt_ideal = 1e10 !Maximum time, will get smaller based on the conditions

    if (proc_num == 0) then
    if (nu0 > 0) dt_ideal = (min(dx,dy,dz))**2/(nu0*(1.0))
    if (nu0 > 0) print*, 'dt due to nu0', cfl*(min(dx,dy,dz))**2/(nu0*(1.0))
    if (voutfact > 0) dt_ideal = min(dt_ideal, dz/voutfact)
    if (voutfact > 0) print*, 'dt due to outflow',  cfl*dz/voutfact
    if (eta > 0 ) dt_ideal = min(dt_ideal, min(dx,dy,dz)**2/eta)
    if (eta > 0 ) print*, 'dt due to eta', cfl*min(dx,dy,dz)**2/eta
    if (shearfact > 0 ) dt_ideal = min(dt_ideal, min(dy,dx)/shearfact)
    if (shearfact > 0 ) print*, 'dt due to shearing', cfl*min(dy,dx)/(shearfact*15.0)
    if (eta0 > 0) dt_ideal = min(dt_ideal, min(dx,dy,dz)**2/eta0)
    if (eta0 > 0) print*, 'dt due to eta0',  cfl*min(dx,dy,dz)**2/eta0
    dt_ideal = cfl*dt_ideal
    !print*, 'Ideal dt', dt_ideal

    dt = dt_ideal

    !print*, 'Final dt', dt, ', total timesteps', nt, ', ', int(nt/(nplots-1)), 'per snapshot'
    end if

    call MPI_Bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    call MPI_Bcast(nt_ft, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)

    nt = int(nt_ft)

    call MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.

    max_velocity = 1D-2*min(dx,dy,dz)/dt

END SUBROUTINE calculate_timestep

SUBROUTINE set_outflow()
    ! Set the outflow arrays vout and voutc, the sizes of which are declared in shared variables
    ! vouts is on gridpoints (like az), voutc is on x ribs (like ax). To allow for upwinding
    ! Need to take care with all this...
    INTEGER:: i, j, k
    REAL(num):: hfact

    allocate(voutx(0:nx+1,-1:ny+1,-1:nz+1))  !Outflow combining with ex
    allocate(vouty(-1:nx+1,0:ny+1,-1:nz+1))  !Outflow combining with ey

    if (voutfact < 0.0) then  !Set to maximum permitted
        voutfact = max_velocity
    end if

    do k = -1, nz+1
        hfact = (zs(k) - z0_global)/(z1_global - z0_global)   !Distance up the domain.
        do i = 0, nx+1
            do j = -1, ny+1
                voutx(i,j,k) = voutfact*hfact**2
            end do
        end do
        do i = -1, nx+1
            do j = 0, ny+1
                vouty(i,j,k) = voutfact*hfact**2
            end do
        end do
    end do

END SUBROUTINE set_outflow

SUBROUTINE set_shearing()
    ! Set the 1D array for the shearing velocity on the lower boundary. Is added to ex so is aligned with the x grid centres
    INTEGER:: i
    allocate(vz0(-1:nx+1))

    do i = -1, nx+1
        vz0(i) = shearfact*sin(xs(i)*PI/x1_global)
    end do

END SUBROUTINE set_shearing

!*******************************************************************************
END MODULE init
!*******************************************************************************
