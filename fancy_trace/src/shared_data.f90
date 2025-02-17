!*******************************************************************************
MODULE shared_data
!*******************************************************************************
! Establish the data that needs to be shared among all the bits of the program
!*******************************************************************************

    IMPLICIT NONE

!*******************************************************************************

    INTEGER, PARAMETER:: d = KIND(0.d0) ! precision for floats
    REAL(d), PARAMETER:: PI = 4.0_d * atan(1.0_d)   ! pi
    INTEGER, PARAMETER :: num = KIND(1.D0)

    ! MPI numbers (total number of processors, etc.)
    INTEGER:: ierr, nprocs, proc_num, comm

    ! Process position information
    REAL(num):: send_test

    !Parameters
    real(num), dimension(:):: flparameters(0:29)
    integer:: run, snap, print_flag, data_source, machine_flag, save_all
    CHARACTER(LEN =128):: data_root
    CHARACTER(LEN =128):: bfield_filename

    !Grid
    integer:: nx, ny, nz
    real(num):: x0, x1, y0, y1, z0, z1
    real(num):: dx, dy, dz

    real(num), dimension(:), allocatable:: xs, ys, zs
    real(num), dimension(:), allocatable:: xc, yc, zc

    !Magnetic field and current
    real(num), dimension(:,:,:), allocatable:: bx, by, bz
    real(num), dimension(:,:,:), allocatable:: jx, jy, jz

    real(num), dimension(:):: b1(0:2), j1(0:2)
    real(num):: jmag

    !Emissivity
    real(num), dimension(:,:,:), allocatable:: emiss
    integer:: nx_out, ny_out, nz_out
    real(num), dimension(:), allocatable:: xs_out, ys_out, zs_out

    !Field line things
    real(num):: ds_factor, t, ds, weakness_limit
    integer:: max_line_length, nstarts, null_point, line_number
    real(num), allocatable:: starts(:,:)

!*******************************************************************************
END MODULE shared_data
!*******************************************************************************
