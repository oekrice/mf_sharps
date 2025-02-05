!*******************************************************************************
MODULE shared_data
!*******************************************************************************
! Establish the data that needs to be shared among all the bits of the program
!*******************************************************************************

    IMPLICIT NONE

    include 'mpif.h'
!*******************************************************************************

    INTEGER, PARAMETER:: d = KIND(0.d0) ! precision for floats
    REAL(d), PARAMETER:: PI = 4.0_d * atan(1.0_d)   ! pi
    INTEGER, PARAMETER :: num = KIND(1.D0)

    !Variable parameters that are read in from the variable text file
    REAL(num):: run_number
    REAL(num):: tmax
    REAL(num):: nplots
    REAL(num):: voutfact
    REAL(num):: bfact
    REAL(num):: shearfact
    REAL(num):: eta
    REAL(num):: eta0
    REAL(num):: nu0
    INTEGER:: init_number
    INTEGER:: mag_min
    INTEGER:: mag_max

    !Other parameters hard-coded into the main.f90 file
    REAL(num):: x0_global, x0
    REAL(num):: x1_global, x1
    REAL(num):: y0_global, y0
    REAL(num):: y1_global, y1
    REAL(num):: z0_global, z0
    REAL(num):: z1_global, z1

    REAL(num):: cfl
    REAL(num):: mf_delta
    CHARACTER(LEN = 64) :: data_directory_root
    CHARACTER(LEN = 64) :: data_directory
    INTEGER:: hamilton_flag
    INTEGER:: ndiags

    !Grid data and arrays shared throughout the code
    INTEGER:: nx_global, nx
    INTEGER:: ny_global, ny
    INTEGER:: nz_global, nz
    INTEGER:: nt, n
    REAL(num):: dx, dy, dz
    REAL(num):: t, dt
    REAL(num), DIMENSION(:), ALLOCATABLE :: xs, ys, zs, xc, yc, zc
    REAL(num), DIMENSION(:), ALLOCATABLE :: xs_global, ys_global, zs_global, xc_global, yc_global, zc_global
    REAL(num):: volume, volume_global

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ax, ay, az
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx, by, bz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jx, jy, jz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ex, ey, ez

    REAL(num), DIMENSION(:,:), ALLOCATABLE :: vx_surf, vy_surf

    REAL(num), DIMENSION(:), ALLOCATABLE:: vz0

    !Temporary arrays that it's probably just quicker to overwrite each time...
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: bx1, by1, bz1
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jx1, jy1, jz1
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: ex1, ey1, ez1

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: vx, vy, vz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: lf, b2, soft, nu

    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: voutx, vouty

    !Extra arrays for the pressure bits
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: fz
    REAL(num):: a, b, deltaz, zstar
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: vpx, vpy, vpz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jpx,jpy
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: jpx1,jpy1
    INTEGER:: decay_type

    !Diagnostics
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_time
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_oflux
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_sumj
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_sume
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_avgj
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_energy
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_maxlorentz
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_avglorentz
    REAL(num), DIMENSION(:), ALLOCATABLE:: diag_nulls
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: diag_lfheights

    !MPI
    INTEGER:: comm, ierr, MPI_loc(3), mpitag = 1
    INTEGER:: nprocs, proc_num
    INTEGER:: x_procs, y_procs, z_procs
    INTEGER:: x_rank, y_rank, z_rank
    INTEGER:: x_up, y_up, z_up, x_down, y_down, z_down
    INTEGER :: bx_xface, by_xface, bz_xface, bx_xface1
    INTEGER :: bx_yface, by_yface, bz_yface, by_yface1
    INTEGER :: bx_zface, by_zface, bz_zface, bz_zface1
    INTEGER :: b0_chunk, b0_all
    INTEGER, DIMENSION(:,:), ALLOCATABLE:: allranks
    !Extras
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: bz_surf_reference
    INTEGER:: nmags
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: surf_vx, surf_vy, surf_vz
    REAL(num), DIMENSION(:,:), ALLOCATABLE:: surf_ex, surf_ey
    REAL(num), DIMENSION(:), ALLOCATABLE:: mag_times

!*******************************************************************************
END MODULE shared_data
!*******************************************************************************
