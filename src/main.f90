program acceleration
    use mpi_f08
    use result
    use user_variables, only: n_sets, filename, outdir, stepsize_exp_str, shockless

    implicit none
    integer :: set
    integer :: myid, n_proc, ierr, n_array
    type(MPI_FILE) :: &
        traj_filehandle, &            ! Used both w/wo shock
        crossang_filehandle, &        ! Only w shock
        phase_space_filehandle, &     ! Only w shock
        fpos_filehandle, &            ! Only wo shock
        samplepos_filehandle          ! Only wo shock
    integer :: traj_array_count, traj_array_bsize
    integer :: crossang_array_count, crossang_array_bsize
    integer :: phase_space_array_count, phase_space_array_bsize
    integer :: fpos_array_count, fpos_array_bsize
    integer :: samplepos_array_count, samplepos_array_bsize
    integer(kind=MPI_OFFSET_KIND) :: &
        traj_offset, &
        crossang_offset, &
        phase_space_offset, &
        fpos_offset, &
        samplepos_offset

    ! MPI init
    call MPI_INIT(ierr)  ! Initiate/create
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) ! Get proc id
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr) ! Get communicator size (# procs)

    ! Simulation init
    call init(myid, n_proc)

    ! MPI data file -- trajectories
    call init_mpi_io(MPI_COMM_WORLD, &
                     traj_filehandle, &
                     trim(outdir)//'/trajectories'//filename, &
                     ierr &
                     )
    if (.not. shockless) then
        ! MPI data file -- crossing flight angles
        call init_mpi_io(MPI_COMM_WORLD, &
                         crossang_filehandle, &
                         trim(outdir)//'/cross_angles'//filename, &
                         ierr &
                         )
        ! MPI data file -- phase space density/distribution
        call init_mpi_io(MPI_COMM_WORLD, &
                         phase_space_filehandle, &
                         trim(outdir)//'/phase_space'//filename, &
                         ierr &
                         )
    else
        ! MPI data file -- final positions
        call init_mpi_io(MPI_COMM_WORLD, &
                         fpos_filehandle, &
                         trim(outdir)//'/fpos'//filename, &
                         ierr &
                         )
        ! MPI data file -- phase space density/distribution
        call init_mpi_io(MPI_COMM_WORLD, &
                         samplepos_filehandle, &
                         trim(outdir)//'/samplepos'//filename, &
                         ierr &
                         )
    end if

    ! Simulation loop -- loop through each set of particles and simulate
    do set = 1, n_sets
        call start_particle(set, myid, n_proc)

        if (.not. shockless) then
            ! Collect data from all procs
            n_array = (2*pid_max + 1)*n_enbin
            call MPI_REDUCE( &
                En_f, &
                En_f_tot, &
                n_array, &
                MPI_DOUBLE_PRECISION, &
                MPI_SUM, &
                0, &
                MPI_COMM_WORLD, &
                ierr &
                )
            call MPI_REDUCE( &
                NE_esc, &
                NE_esc_tot, &
                n_enbin, &
                MPI_DOUBLE_PRECISION, &
                MPI_SUM, &
                0, &
                MPI_COMM_WORLD, &
                ierr &
                )
            if (myid == 0) call output(set, n_proc)
        end if

        ! Write trajectory data
        traj_array_bsize = sizeof(trajectories)
        traj_offset = myid*n_sets*traj_array_bsize + traj_array_bsize*(set - 1)
        traj_array_count = size(trajectories)
        call MPI_FILE_WRITE_AT( &
            traj_filehandle, &
            traj_offset, &
            trajectories, &
            traj_array_count, &
            MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE, &
            ierr &
            )
        if (.not. shockless) then
            ! Write cross angle data
            crossang_array_bsize = sizeof(crossing_flight_angles)
            crossang_offset = myid*n_sets*traj_array_bsize + traj_array_bsize*(set - 1)
            crossang_array_count = size(crossing_flight_angles)
            call MPI_FILE_WRITE_AT( &
                crossang_filehandle, &
                crossang_offset, &
                crossing_flight_angles, &
                crossang_array_count, &
                MPI_DOUBLE_PRECISION, &
                MPI_STATUS_IGNORE, &
                ierr &
                )
            ! Write phase data
            phase_space_array_bsize = sizeof(phase_space_pos)
            phase_space_offset = myid*n_sets*phase_space_array_bsize + phase_space_array_bsize*(set - 1)
            phase_space_array_count = size(phase_space_pos)
            call MPI_FILE_WRITE_AT( &
                phase_space_filehandle, &
                phase_space_offset, &
                phase_space_pos, &
                phase_space_array_count, &
                MPI_DOUBLE_PRECISION, &
                MPI_STATUS_IGNORE, &
                ierr &
                )
        else
            ! Write fpos data
            fpos_array_bsize = sizeof(final_positions)
            fpos_offset = myid*n_sets*fpos_array_bsize + fpos_array_bsize*(set - 1)
            fpos_array_count = size(final_positions)
            call MPI_FILE_WRITE_AT( &
                fpos_filehandle, &
                fpos_offset, &
                final_positions, &
                fpos_array_count, &
                MPI_DOUBLE_PRECISION, &
                MPI_STATUS_IGNORE, &
                ierr &
                )
            ! Write sample_positions
            samplepos_array_bsize = sizeof(sample_positions)
            samplepos_offset = &
                myid*n_sets*samplepos_array_bsize + samplepos_array_bsize*(set - 1)
            samplepos_array_count = size(sample_positions)
            call MPI_FILE_WRITE_AT( &
                samplepos_filehandle, &
                samplepos_offset, &
                sample_positions, &
                samplepos_array_count, &
                MPI_DOUBLE_PRECISION, &
                MPI_STATUS_IGNORE, &
                ierr &
                )
        end if
    end do

    ! close outputs and finalize
    call MPI_FILE_CLOSE(traj_filehandle, ierr)
    if (.not. shockless) then
        call MPI_FILE_CLOSE(crossang_filehandle, ierr)
        call MPI_FILE_CLOSE(phase_space_filehandle, ierr)
    else
        call MPI_FILE_CLOSE(fpos_filehandle, ierr)
        call MPI_FILE_CLOSE(samplepos_filehandle, ierr)
    end if
    close (99)
    call MPI_FINALIZE(ierr)
end program acceleration

subroutine start_particle(set, myid, n_proc)
    use user_variables, only: n_start!, debug
    use internal, only: n_in
    use test_var, only: n_injected, sec

    implicit none
    integer set, myid, n_proc

    n_injected = 0
    do while (n_injected < n_start)
        if (n_in == 0) then
            n_injected = n_injected + 1
            call inject !(n_injected)
            n_in = 1
            sec = 0
            ! write(*,*)
            ! write(*,*) 'primary',n_injected,n_start
        else
            sec = 1
            ! write(*,*) 'secondary'
        end if
        call tracer(set, n_injected)
        if (myid == 0 .and. mod(n_injected*100, n_start) == 0 .and. sec == 0) &
            write (*, *) set, n_injected*n_proc
    end do
end subroutine start_particle

subroutine tracer(set, n_injected)
    use event_internal; use internal, only: n_in
    use random_walk
    use user_variables, only: isotropic

    implicit none
    integer id
    integer, pointer :: A, pid
    integer, intent(in) :: set, n_injected

    A => event(n_in)%A
    pid => event(n_in)%pid
    if (A > 1) then
        id = 100 + A
    else
        id = abs(pid)
    end if
    select case (id)
        ! case (102:108)                                      ! discard low A nuclei
    case (102:144)                                        ! discard low A nuclei
        n_in = n_in - 1
        return
    case (7, 145:159)
        ! Start particle simulation
        call pitch_angle_random_walk(set, n_injected)
    case default
        write (*, *) 'A,pid', A, pid
        call error('wrong particle typ in tracer', 0)
    end select
end subroutine tracer
