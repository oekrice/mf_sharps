!*******************************************************************************
PROGRAM main
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code
!*******************************************************************************
    USE shared_data
    USE init
    USE mpi_tools
    USE evolve
    USE output
    !USE output, ONLY: save_snap, print_array, diagnostics
    IMPLICIT NONE

    INTEGER:: nstart, nend, init_mag, mag_interval
    INTEGER:: block_num, diag_num
    REAL(num):: tmags, block_dt, tstart, tend, diag_ideal_time
    LOGICAL:: first_diagnostic
    !REAL(num):: mag_interval
    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs
    cfl  = 0.1
    mf_delta = 1D-3

    ! Import the parameters and set up the grid
    CALL initialise()

    print*, 'timestep', dt

    if (.true.) then
    if (hamilton_flag < 0.5) then
        data_directory_root = '/extra/tmp/trcn27/mf3d/'
    else
        data_directory_root = '/nobackup/trcn27/mf3d0/'
    end if

    if (run_number < 10) then
        write (data_directory, "(A23,A2,I1,A1)") data_directory_root, '00',int(run_number), "/"
    else if (run_number < 100) then
        write (data_directory, "(A23,A1,I2,A1)") data_directory_root, '0', int(run_number), "/"
    else
        write (data_directory, "(A23,I3,A1)") data_directory_root, int(run_number), "/"
    end if

    !CALL update_surface_flows(0)

    if (proc_num == 0) print*, 'Set up, running from', mag_min, 'to', mag_max
    if (proc_num == 0) print*, 'Times from', mag_times(mag_min), 'to', mag_times(mag_max)

    if (mag_min == 0) then
        CALL timestep()  !Does everything except the actual timestep (for diagnostic reasons)
        CALL export_magnetogram(0)
        CALL save_snap(0)
    end if

    !Calculate initial snapshot and diagnostic numbers -- don't have to line up with magnetograms
    first_diagnostic = .true.
    do block_num = mag_min, mag_max-1
        tstart = mag_times(block_num) ; tend = mag_times(block_num+1)
        CALL import_surface_electric(block_num, 1.0_num/(tend - tstart))
        CALL import_surface_magnetic(block_num)
        t = tstart
        block_dt = (tend-tstart)/(int((tend-tstart)/dt) + 1)
        nt = int((tend-tstart)/block_dt)
        !Within each block, readjust timesteps and things
        if (tstart < 1D-6) then
            diag_num = 0
        else
            diag_num = int(ndiags*tstart/tmax) + 1  !Initial diagnostic number
        end if
        diag_ideal_time = diag_num*tmax/ndiags

        do n = 0, nt-1  ! Actually run the code
            mag_ratio = (float(n) + 0.5_num)/nt
            CALL timestep()  !Does everything except the actual timestep (for diagnostic reasons)
            !Check whether diagnostics or a snapshot is necessary
            if (t - block_dt < diag_ideal_time .and. t .ge. diag_ideal_time) then
                CALL diagnostics(diag_num, first_diagnostic)
                diag_num = diag_num + 1
                diag_ideal_time = diag_num*tmax/ndiags
                first_diagnostic = .false.
            end if
            ax = ax - block_dt*ex
            ay = ay - block_dt*ey
            az = az - block_dt*ez

            t = t + block_dt

        end do

        !Need an extra 'timestep' here to get the magnetic field right for the export -- but don't update A or the time
        CALL timestep()
        CALL export_magnetogram(block_num+1)
        CALL save_snap(block_num+1)

        if (t - block_dt < diag_ideal_time .and. t .ge. diag_ideal_time) then
            CALL diagnostics(diag_num, first_diagnostic)
            diag_num = diag_num + 1
            diag_ideal_time = diag_num*tmax/ndiags
            first_diagnostic = .false.
        end if

        if (proc_num == 0) print*, 'Snap saved', block_num+1, t

    end do

    if (proc_num == -1) then
        print*, 'Open Flux', proc_num, z_rank, sum(abs(bz(1:nx,1:ny,nz)))
        print*, 'Max. currents', proc_num, sum(abs(jx(2:nx-2,2:ny-2,2:nz-1))), &
        sum(abs(jy(2:nx-2,2:ny-2,2:nz-1))), sum(abs(jz(2:nx-2,2:ny-2,2:nz-1)))
    end if
    !CALL diagnostics(int(n/(nt/(ndiags-1))))
    if (proc_num == -1) print*, 'Step', n, 'at time', t

    !CALL diagnostics(ndiags-1)
    end if
    !if (proc_num == 0 .and. mag_max == 500) print*, 'Fortran code completed sucessfully. Carry on.'
    CALL mpi_finalize(ierr)
    STOP

END PROGRAM main
