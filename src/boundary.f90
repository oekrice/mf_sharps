!*******************************************************************************

MODULE boundary
!*******************************************************************************
! Initialise the simulation. Read in/compute the initial condition (python code for that I think?), read in the parameters and all will be well
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    !USE output
    !USE pressure
    IMPLICIT NONE

!*******************************************************************************
CONTAINS

SUBROUTINE magnetic_boundary()

    IMPLICIT NONE

    CALL bfield_mpi()

    !If first step, establish reference magnetic field

    CALL MPI_BARRIER(comm,ierr)

    !Boundary conditions on the magnetic field.
    !Uses the no-horizontal-current condition for the first ghost cells and then the solenoidal condition for the next lot.

    !PROBLEM WITH CORNERS!


    !LOWER BOUNDARY (Zero current)
    if (z_rank == 0) then
    surf_by(0:nx+1,0:ny) = mag_ratio*surf_by1(0:nx+1,0:ny) + (1.0_num - mag_ratio)*surf_by0(0:nx+1,0:ny)
    surf_bx(0:nx, 0:ny+1) = mag_ratio*surf_bx1(0:nx, 0:ny+1) + (1.0_num - mag_ratio)*surf_bx0(0:nx, 0:ny+1)

    by(0:nx+1,0:ny,0) = surf_by(0:nx+1,0:ny)
    bx(0:nx,0:ny+1,0) = surf_bx(0:nx,0:ny+1)

!     by(0:nx+1,0:ny,0) = by(0:nx+1,0:ny,1) - dz*(bz(0:nx+1,1:ny+1,0) - bz(0:nx+1, 0:ny,0))/dy
!     bx(0:nx, 0:ny+1,0) = bx(0:nx,0:ny+1,1) - dz*(bz(1:nx+1,0:ny+1,0) - bz(0:nx,0:ny+1,0))/dx
    end if

    !UPPER BOUNDARY (Zero Current)
    if (z_rank == z_procs-1) then
    by(0:nx+1,0:ny,nz+1) = by(0:nx+1,0:ny,nz) + dz*(bz(0:nx+1,1:ny+1,nz) - bz(0:nx+1, 0:ny,nz))/dy
    bx(0:nx, 0:ny+1,nz+1) = bx(0:nx,0:ny+1,nz) + dz*(bz(1:nx+1,0:ny+1,nz) - bz(0:nx,0:ny+1,nz))/dx
    end if

      !x boundaries (Zero current, and zero flux)
    if (x_rank == 0) then
      bx(-1,:,:) = bx(1,:,:)
      by( 0,:,:) = -by(1,:,:)
      bz( 0,:,:) = bz(1,:,:)
    end if

    if (x_rank == x_procs-1) then
      bx(nx+1,:,:) = bx(nx-1,:,:)
      by(nx+1,:,:) = -by(nx  ,:,:)
      bz(nx+1,:,:) = bz(nx  ,:,:)
    end if

    !y boundaries (Zero current, and zero flux)
    if (y_rank == 0) then
      bx(:, 0,:) = -bx(:,1,:)
      by(:,-1,:) = by(:,1,:)
      bz(:, 0,:) = bz(:,1,:)
    end if

    if (y_rank == y_procs-1) then

      bx(:,ny+1,:) = -bx(:,ny  ,:)
      by(:,ny+1,:) = by(:,ny-1,:)
      bz(:,ny+1,:) = bz(:,ny  ,:)

    end if

    if (z_rank == 0) then
      bx(:,:, 0) = bx(:,:,1)
      by(:,:, 0) = by(:,:,1)
      bz(:,:,-1) = bz(:,:,1)

    end if

    !UPPER BOUNDARY (Zero Current)
    if (z_rank == z_procs-1) then
      bx(:,:,nz+1) = bx(:,:,nz  )
      by(:,:,nz+1) = by(:,:,nz  )
      bz(:,:,nz+1) = bz(:,:,nz-1)
    end if

    !x boundaries (Zero current, and zero flux)
    if (x_rank == 0) then
    bz(0,0:ny+1,0:nz) = bz(1,0:ny+1,0:nz) - dx*(bx(0,0:ny+1,1:nz+1) - bx(0, 0:ny+1,0:nz))/dz
    !by(0,0:ny,0:nz+1) = by(1, 0:ny,0:nz+1) - dx*(bx(0,1:ny+1,0:nz+1) - bx(0,0:ny,0:nz+1))/dy
    end if

    if (x_rank == x_procs-1) then
    bz(nx+1,0:ny+1,0:nz) = bz(nx,0:ny+1,0:nz) + dx*(bx(nx,0:ny+1,1:nz+1) - bx(nx, 0:ny+1,0:nz))/dz
    !by(nx+1,0:ny,0:nz+1) = by(nx, 0:ny,0:nz+1) + dx*(bx(nx,1:ny+1,0:nz+1) - bx(nx,0:ny,0:nz+1))/dy
    end if

    !y boundaries (Zero current, and zero flux)
    if (y_rank == 0) then
    bz(0:nx+1,0,0:nz) = bz(0:nx+1, 1,0:nz) - dy*(by(0:nx+1,0,1:nz+1) - by(0:nx+1,0,0:nz))/dz
    !bx(0:nx,0,0:nz+1) = bx(0:nx,1,0:nz+1) - dy*(by(1:nx+1,0,0:nz+1) - by(0:nx, 0,0:nz+1))/dx
    end if

    if (y_rank == y_procs-1) then
    bz(0:nx+1,ny+1,0:nz) = bz(0:nx+1, ny,0:nz) + dy*(by(0:nx+1,ny,1:nz+1) - by(0:nx+1,ny,0:nz))/dz
    !bx(0:nx,ny+1,0:nz+1) = bx(0:nx,ny,0:nz+1) + dy*(by(1:nx+1,ny,0:nz+1) - by(0:nx, ny,0:nz+1))/dx
    end if

    CALL MPI_BARRIER(comm,ierr)

    !_________________________________________________
    !Solenoidal condition, using the above information
    !Isn't necessary to do more MPI here, I suppose. So I won't.

    !LOWER BOUNDARY
    if (z_rank == 0) then
    bz(0:nx+1,0:ny+1,-1) = (1.0_num/(dx*dy))*(bz(0:nx+1,0:ny+1,0)*dx*dy - bx(-1:nx,0:ny+1,0)*dy*dz + &
    bx(0:nx+1,0:ny+1,0)*dy*dz - by(0:nx+1,-1:ny,0)*dx*dz + by(0:nx+1,0:ny+1,0)*dx*dz)
    end if

    !UPPER BOUNDARY
    if (z_rank == z_procs - 1) then
    bz(0:nx+1,0:ny+1,nz+1) = (1.0_num/(dx*dy))*(bz(0:nx+1,0:ny+1,nz)*dx*dy + bx(-1:nx,0:ny+1,nz+1)*dy*dz - &
    bx(0:nx+1,0:ny+1,nz+1)*dy*dz + by(0:nx+1,-1:ny,nz+1)*dx*dz - by(0:nx+1,0:ny+1,nz+1)*dx*dz)
    end if

    !x boundaries (Zero current, and zero flux)
    if (x_rank == 0) then
    bx(-1,0:ny+1,0:nz+1) = (1.0_num/(dy*dz))*(bx(0,0:ny+1,0:nz+1)*dy*dz - by(0,-1:ny,0:nz+1)*dx*dz + &
    by(0,0:ny+1,0:nz+1)*dx*dz - bz(0,0:ny+1,-1:nz)*dx*dy + bz(0,0:ny+1,0:nz+1)*dx*dy)
    end if

    if (x_rank == x_procs-1) then
    bx(nx+1,0:ny+1,0:nz+1) = (1.0_num/(dy*dz))*(bx(nx,0:ny+1,0:nz+1)*dy*dz + by(nx+1,-1:ny,0:nz+1)*dx*dz - &
    by(nx+1,0:ny+1,0:nz+1)*dx*dz + bz(nx+1,0:ny+1,-1:nz)*dx*dy - bz(nx+1,0:ny+1,0:nz+1)*dx*dy)
    end if

    !y boundaries (Zero current, and zero flux)
    if (y_rank == 0) then
    by(0:nx+1,-1,0:nz+1) = (1.0_num/(dx*dz))*(by(0:nx+1,0,0:nz+1)*dx*dz - bx(-1:nx,0,0:nz+1)*dy*dz + &
    bx(0:nx+1,0,0:nz+1)*dy*dz - bz(0:nx+1,0,-1:nz)*dx*dy + bz(0:nx+1,0,0:nz+1)*dx*dy)
    end if

    if (y_rank == y_procs-1) then
    by(0:nx+1,ny+1,0:nz+1) = (1.0_num/(dx*dz))*(by(0:nx+1,ny,0:nz+1)*dx*dz + bx(-1:nx,ny+1,0:nz+1)*dy*dz - &
    bx(0:nx+1,ny+1,0:nz+1)*dy*dz + bz(0:nx+1,ny+1,-1:nz)*dx*dy - bz(0:nx+1,ny+1,0:nz+1)*dx*dy)
    end if

END SUBROUTINE magnetic_boundary

!*******************************************************************************
END MODULE boundary
!*******************************************************************************
