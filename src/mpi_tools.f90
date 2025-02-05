!*******************************************************************************
MODULE mpi_tools
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code. Theoretically. Want to run with intel debuggers on, otherwise definitely not going to pick up on all the mistakes.
!*******************************************************************************
    USE shared_data

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE start_mpi()

        IMPLICIT NONE
        integer:: i
        INTEGER:: mpi_dims(3), mpi_dims_best(3)
        LOGICAL:: mpi_periodic(3)
        REAL(num):: diff_best, diff, mean

        ! - Have already established the global rank.
        call mpi_init(ierr)  !Tells it to start using MPI

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) !Number of processes globally.
        call mpi_comm_rank(MPI_COMM_WORLD, proc_num, ierr) !Returns the rank of current process

        ! Choose optimum division of procs that fits grid dimensions:
        ! This algorithm only appears to work for 4 or more processors, so will need some exceptions
        if (nprocs == 1) then
            mpi_dims_best = (/1,1,1/)

        else if (nprocs == 2) then
            mpi_dims_best = (/1,1,2/)
        else
        diff_best = REAL(nprocs)
        mpi_dims_best = (/-1,-1,-1/)
        DO i = 2, nprocs, 2
            IF (MOD(nprocs, i) == 0) THEN
                mpi_dims = (/0, 0, i/)
                ! Find optimum decomposition with i points in p:
                CALL MPI_DIMS_CREATE(nprocs, 3, mpi_dims, ierr)
                ! Check whether this is allowed:
                nx = nx_global / mpi_dims(1)
                ny = ny_global / mpi_dims(2)
                nz = nz_global / mpi_dims(3)
                IF (((nx * mpi_dims(1)) == nx_global) &
                    .AND. ((ny * mpi_dims(2) == ny_global)) &
                    .AND. ((nz * mpi_dims(3)) == nz_global)) THEN
                    mean = SUM(REAL(mpi_dims))/3.0_d
                    diff = MAXVAL(REAL(mpi_dims) - mean)
                    IF (diff < diff_best) THEN
                        diff_best = diff
                        mpi_dims_best = mpi_dims
                    END IF
                END IF
            END IF
        END DO
        end if
        IF (mpi_dims_best(1) * mpi_dims_best(2) * mpi_dims_best(3) == nprocs) THEN
            mpi_dims = mpi_dims_best
            nx = nx_global / mpi_dims(1)
            ny = ny_global / mpi_dims(2)
            nz = nz_global / mpi_dims(3)
        ELSE
            PRINT*,'ERROR: THIS NUMBER OF MPI PROCS DOES NOT FIT THE GRID'
            CALL MPI_abort(comm, ierr)
        END IF

        x_procs = mpi_dims(1); y_procs = mpi_dims(2); z_procs = mpi_dims(3);
        MPI_periodic = (/.false.,.false.,.false./)
        !Attempt to use the DUMFRIC way of establishing the communicator:
        CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, MPI_dims, MPI_periodic, .TRUE., &
        comm, ierr)

        !Redistribute ranks based on this new communicator
        CALL MPI_COMM_RANK(comm, proc_num, ierr)

        CALL MPI_CART_COORDS(comm, proc_num, 3, mpi_loc, ierr)
        CALL MPI_CART_SHIFT(comm, 0, 1, x_down, x_up, ierr)
        CALL MPI_CART_SHIFT(comm, 1, 1, y_down, y_up, ierr)
        CALL MPI_CART_SHIFT(comm, 2, 1, z_down, z_up, ierr)

        x_rank = mpi_loc(1); y_rank = mpi_loc(2); z_rank = mpi_loc(3)

        call MPI_BARRIER(comm, ierr)

        do i = 0, nprocs-1
            if (proc_num == i .and. .false.) then
                print*, proc_num, x_rank, y_rank, z_rank
                !print*, 'x', x_down, x_up
                !print*, 'y', y_down, y_up
                !print*, 'z', z_down, z_up
                print*, '______________________________'
            end if
            call MPI_BARRIER(comm, ierr)
        end do

        !Distribute MPI rank information, for snazzy diagnostics
        allocate(allranks(0:nprocs-1, 0:2))

        do i = 0, nprocs-1
        if (i == proc_num) then
        call MPI_ALLREDUCE((/x_rank, y_rank, z_rank/),allranks(i,0:2),3,MPI_INTEGER,MPI_SUM,comm,ierr)
        else
        call MPI_ALLREDUCE((/0, 0, 0/),allranks(i,0:2),3,MPI_INTEGER,MPI_SUM,comm,ierr)
        end if
        end do

        return
    END SUBROUTINE start_mpi

    SUBROUTINE bfield_mpi

        IMPLICIT NONE

        !MPI routines for passing the boundary data
        !Send z data DOWN
        if (z_down >= 0) then
            call mpi_send(by(0,-1,1), 1, by_zface, z_down, 0, comm, ierr)
            call mpi_send(bx(-1,0,1), 1, bx_zface, z_down, 1, comm, ierr)
        end if

        !Receive z data from ABOVE
        if (z_up >= 0) then
            call mpi_recv(by(0,-1,nz+1), 1, by_zface, z_up, 0, comm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(bx(-1,0,nz+1), 1, bx_zface, z_up, 1, comm, MPI_STATUS_IGNORE, ierr)
        end if

        !Send z data UP
        if (z_up >= 0) then
            call mpi_send(by(0,-1,nz), 1, by_zface, z_up, 0, comm, ierr)
            call mpi_send(bx(-1,0,nz), 1, bx_zface, z_up, 1, comm, ierr)
        end if

        !Receive x data from BELOW
        if (z_down >= 0) then
            call mpi_recv(by(0,-1,0), 1, by_zface, z_down, 0, comm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(bx(-1,0,0), 1, bx_zface, z_down, 1, comm, MPI_STATUS_IGNORE, ierr)
        end if

        if (x_down >= 0) then
            call mpi_send(by(1,-1,0), 1, by_xface, x_down, 0, comm, ierr)
            call mpi_send(bz(1,0,-1), 1, bz_xface, x_down, 1, comm, ierr)
        end if

        !Receive x data from ABOVE
        if (x_up >= 0) then
            call mpi_recv(by(nx+1,-1,0), 1, by_xface, x_up, 0, comm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(bz(nx+1,0,-1), 1, bz_xface, x_up, 1, comm, MPI_STATUS_IGNORE, ierr)
        end if

        !Send x data UP
        if (x_up >= 0) then
            call mpi_send(by(nx,-1,0), 1, by_xface, x_up, 0, comm, ierr)
            call mpi_send(bz(nx,0,-1), 1, bz_xface, x_up, 1, comm, ierr)
        end if

        !Receive x data from BELOW
        if (x_down >= 0) then
            call mpi_recv(by(0,-1,0), 1, by_xface, x_down, 0, comm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(bz(0,0,-1), 1, bz_xface, x_down, 1, comm, MPI_STATUS_IGNORE, ierr)
        end if

        if (y_down >= 0) then
            call mpi_send(bx(-1,1,0), 1, bx_yface, y_down, 0, comm, ierr)
            call mpi_send(bz(0,1,-1), 1, bz_yface, y_down, 1, comm, ierr)
        end if

        !Receive y data from ABOVE
        if (y_up >= 0) then
            call mpi_recv(bx(-1,ny+1,0), 1, bx_yface, y_up, 0, comm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(bz(0,ny+1,-1), 1, bz_yface, y_up, 1, comm, MPI_STATUS_IGNORE, ierr)
        end if

        !Send y data UP
        if (y_up >= 0) then
            call mpi_send(bx(-1,ny,0), 1, bx_yface, y_up, 0, comm, ierr)
            call mpi_send(bz(0,ny,-1), 1, bz_yface, y_up, 1, comm, ierr)
        end if

        !Receive y data from BELOW
        if (y_down >= 0) then
            call mpi_recv(bx(-1,0,0), 1, bx_yface, y_down, 0, comm, MPI_STATUS_IGNORE, ierr)
            call mpi_recv(bz(0,0,-1), 1, bz_yface, y_down, 1, comm, MPI_STATUS_IGNORE, ierr)
        end if


        call MPI_BARRIER(comm, ierr)

        return
    END SUBROUTINE bfield_mpi

    SUBROUTINE mpi_create_types
        !There appears to be a (fairly low) integer limit on the MPI, so I'll try to get around this by doing the types properly.
        !Should hopefully only need to do six of these.
        IMPLICIT NONE
        INTEGER:: mpitype

        !Types for the sends in the magnetic boundary conditions
        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx+2,ny+3,nz+2/), (/nx+2,ny+3,1/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        by_zface = mpitype

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx+3,ny+2,nz+2/), (/nx+3,ny+2,1/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        bx_zface = mpitype

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx+3,ny+2,nz+2/), (/nx+3,1,nz+2/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        bx_yface = mpitype

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx+2,ny+2,nz+3/), (/nx+2,1,nz+3/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        bz_yface = mpitype

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx+2,ny+3,nz+2/), (/1,ny+3,nz+2/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        by_xface = mpitype

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx+2,ny+2,nz+3/), (/1,ny+2,nz+3/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        bz_xface = mpitype

        !Types for sending an entire chunk of averaged field

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx_global,ny_global,nz_global/), (/nx,ny,nz/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        b0_chunk = mpitype

        !Types for sending an entire averaged field

        mpitype = MPI_DATATYPE_NULL
        CALL MPI_TYPE_CREATE_SUBARRAY(3, (/nx_global, ny_global, nz_global/), (/nx_global,ny_global,nz_global/), (/0,0,0/), &
            MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ierr)
        CALL MPI_TYPE_COMMIT(mpitype, ierr)

        b0_all = mpitype



    END SUBROUTINE mpi_create_types

END MODULE mpi_tools

