!*******************************************************************************
MODULE grid
!*******************************************************************************
! Generates all the shared grid data, as grid3d does in the python code. Maybe additionally put shared grid in arrays in a different file - not sure yet. Most (if not all) of the grid arrays should not depend on the process, but I guess we'll see about that.
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

    contains

    subroutine establish_grid()
        !Establishes the grid arrays, which hopefully can be read in mainly from the netcdf. Perhaps.
        !Also read in the magnetic field here because may as well
        IMPLICIT NONE
        INTEGER:: ncid, vid
        character(len=4) :: snap_id
        call get_coordinates

        allocate(bx(0:nx,-1:ny+1,-1:nz+1))
        allocate(by(-1:nx+1,0:ny,-1:nz+1))
        allocate(bz(-1:nx+1,-1:ny+1,0:nz))

        allocate(jx(0:nx+1,0:ny  ,0:nz  ))
        allocate(jy(0:nx  ,0:ny+1,0:nz  ))
        allocate(jz(0:nx  ,0:ny  ,0:nz+1))

        allocate(emiss(0:nx_out,0:ny_out,0:nz_out))

        emiss = 0.0_num

        write (snap_id,'(I4.4)') snap
        bfield_filename = trim(trim(data_root)//trim(snap_id)//'.nc')

        if (print_flag > 0.5) print*, 'Using magnetic field from file ', bfield_filename

        call try(nf90_open(trim(bfield_filename), nf90_nowrite, ncid))

        call try(nf90_inq_varid(ncid, 'bx', vid))
        call try(nf90_get_var(ncid, vid, bx(0:nx,1:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'by', vid))
        call try(nf90_get_var(ncid, vid, by(1:nx,0:ny,1:nz)))

        call try(nf90_inq_varid(ncid, 'bz', vid))
        call try(nf90_get_var(ncid, vid, bz(1:nx,1:ny,0:nz)))

        call try(nf90_inq_varid(ncid, 'jx', vid))
        call try(nf90_get_var(ncid, vid, jx(1:nx,0:ny,0:nz)))

        call try(nf90_inq_varid(ncid, 'jy', vid))
        call try(nf90_get_var(ncid, vid, jy(0:nx,1:ny,0:nz)))

        call try(nf90_inq_varid(ncid, 'jz', vid))
        call try(nf90_get_var(ncid, vid, jz(0:nx,0:ny,1:nz)))

        call try(nf90_close(ncid))

        !call check_divergence

        if (print_flag > 0.5) print*, 'Grid established and magnetic field read-in'

        ds = ds_factor*min(dx, dy, dz)  !Tracing 'timestep'

    end subroutine establish_grid

    subroutine get_coordinates
    !Allocate dimension arrays
    IMPLICIT NONE
    INTEGER:: i,j,k
    allocate(xs(0:nx),ys(0:ny),zs(0:nz))
    allocate(xc(1:nx),yc(1:ny),zc(1:nz))
    allocate(xs_out(0:nx_out),ys_out(0:ny_out),zs_out(0:nz_out))

    dx = (x1 - x0)/nx;     dy = (y1 - y0)/ny;     dz = (z1 - z0)/nz

    xs(0) = x0; ys(0) = y0; zs(0) = z0
    do i = 1, nx
        xs(i) = xs(i-1) + dx
    end do
    do j = 1, ny
        ys(j) = ys(j-1) + dy
    end do
    do k = 1, nz
        zs(k) = zs(k-1) + dz
    end do

    xc(1:nx) = 0.5_num*(xs(0:nx-1) + xs(1:nx))
    yc(1:ny) = 0.5_num*(ys(0:ny-1) + ys(1:ny))
    zc(1:nz) = 0.5_num*(zs(0:nz-1) + zs(1:nz))

    dx = sum((xs(1:nx) - xs(0:nx-1)))/ nx
    dy = sum((ys(1:ny) - ys(0:ny-1)))/ ny
    dz = sum((zs(1:nz) - zs(0:nz-1)))/ nz

    do i = 0, nx_out
        xs_out(i) = xs(0) + (xs(nx) - xs(0))*i/nx_out
    end do

    do j = 0, ny_out
        ys_out(j) = ys(0) + (ys(ny) - ys(0))*j/ny_out
    end do

    do k = 0, nz_out
        zs_out(k) = zs(0) + (zs(nz) - zs(0))*k/nz_out
    end do

    end subroutine get_coordinates

    subroutine check_divergence
    !Checks that the solenoidal condition is fine in every cell
    IMPLICIT none
    real(num):: div(1:nx,1:ny,1:nz)

    div = 0.0_num
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (bx(1:nx,1:ny,1:nz) - bx(0:nx-1,1:ny,1:nz))/dx
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (by(1:nx,1:ny,1:nz) - by(1:nx,0:ny-1,1:nz))/dy
    div(1:nx,1:ny,1:nz) = div(1:nx,1:ny,1:nz) + (bz(1:nx,1:ny,1:nz) - bz(1:nx,1:ny,0:nz-1))/dz

    if (maxval(abs(div)) > 1d-10) then
        print*, 'Read-in magnetic field is not divergence free, stopping'
        STOP
    end if

    end subroutine check_divergence

    SUBROUTINE export_emissivity
    IMPLICIT NONE

    character(len=64):: filename
    integer:: aid, bid, cid, vid, ncid
    integer:: xid, yid, zid
    CHARACTER(LEN=4):: snap_id

    real(num), dimension(:,:):: xsum(0:ny_out,0:nz_out), ysum(0:nx_out,0:nz_out), zsum(0:nx_out,0:ny_out)

    write (snap_id,'(I4.4)') snap

    filename = trim('./fl_data/emiss'//trim(snap_id)//'.nc')

    call try(nf90_create(trim(filename), nf90_clobber, ncid))

    !Define variables
    call try(nf90_def_dim(ncid, 'nx_out', nx_out+1, aid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'ny_out', ny_out+1, bid))  !Make up fake dimensions here
    call try(nf90_def_dim(ncid, 'nz_out', nz_out+1, cid))  !Make up fake dimensions here

    if (save_all > 0.5_num) call try(nf90_def_var(ncid, 'emiss', nf90_double, (/aid, bid, cid/), vid))
    call try(nf90_def_var(ncid, 'emiss_xsum', nf90_double, (/bid, cid/), xid))
    call try(nf90_def_var(ncid, 'emiss_ysum', nf90_double, (/aid, cid/), yid))
    call try(nf90_def_var(ncid, 'emiss_zsum', nf90_double, (/aid, bid/), zid))

    call try(nf90_enddef(ncid))

    xsum = sum(emiss, dim = 1)
    ysum = sum(emiss, dim = 2)
    zsum = sum(emiss, dim = 3)

    !Write variables
    if (save_all > 0.5_num) call try(nf90_put_var(ncid, vid, emiss))
    call try(nf90_put_var(ncid, xid, xsum))
    call try(nf90_put_var(ncid, yid, ysum))
    call try(nf90_put_var(ncid, zid, zsum))

    call try(nf90_close(ncid))

    if (print_flag > 0.5_num) print*, 'Emissivity exported to file', filename

    END SUBROUTINE export_emissivity


    SUBROUTINE try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        print*, 'Ensure data directory is correct in fltrace.f90'
        stop
    end if

    END SUBROUTINE try

END MODULE grid
