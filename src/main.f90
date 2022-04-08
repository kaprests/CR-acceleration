! File: main.f90
! Cosmic ray acceleration

program acceleration
    use mpi_f08
    use user_variables
    use result; use user_variables, only: n_sets

    implicit none
    integer :: myid, n_proc, ierr, n_array
    integer :: set
    type(MPI_FILE) ::  &
        traj_filehandle, &
        fpos_filehandle, &
        samplepos_filehandle
    integer :: traj_array_count, traj_array_bsize
    integer :: fpos_array_count, fpos_array_bsize
    integer :: samplepos_array_count, samplepos_array_bsize
    integer(kind=MPI_OFFSET_KIND) :: &
        traj_offset, &
        fpos_offset, &
        samplepos_offset
    ! non-MPI values
    !myid = 0
    !n_proc = 1

    ! MPI init
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)

    ! Initialize simulation
    call init(myid, n_proc)

    ! MPI data files
    call init_mpi_io(MPI_COMM_WORLD, traj_filehandle, &
        trim(outdir)//'/trajectories'//filename, ierr)
    if (shockless) then
        call init_mpi_io(MPI_COMM_WORLD, fpos_filehandle, &
            trim(outdir)//'/fpos'//filename, ierr)
        call init_mpi_io(MPI_COMM_WORLD, samplepos_filehandle, &
            trim(outdir)//'/samplepos'//filename, ierr)
    end if

    do set = 1, n_sets
        call start_particle(set, myid, n_proc)

        ! non-MPI values
        !En_f_tot = En_f

        ! Write trajectory data
        traj_array_bsize = sizeof(initial_trajectories)
        traj_offset = myid*n_sets*traj_array_bsize + traj_array_bsize*(set - 1)
        traj_array_count = size(initial_trajectories)
        call MPI_FILE_WRITE_AT( &
            traj_filehandle, &
            traj_offset, &
            initial_trajectories, &
            traj_array_count, &
            MPI_DOUBLE_PRECISION, &
            MPI_STATUS_IGNORE, &
            ierr &
            )
        if (shockless) then
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
        else
            ! Collect data from the procs
            n_array = (2*pid_max + 1)*n_enbin
            call MPI_REDUCE( &
                En_f, En_f_tot, n_array, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                MPI_COMM_WORLD, ierr &
                )
            call MPI_REDUCE( &
                cross_angle_distribution, cross_angle_distribution_tot, &
                n_angle_bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
                MPI_COMM_WORLD, ierr &
                )
            if (myid == 0) call output(set, n_proc)
        end if
    end do
    ! close outputs and finalize
    call MPI_FILE_CLOSE(traj_filehandle, ierr)
    if (shockless) then
        call MPI_FILE_CLOSE(fpos_filehandle, ierr)
        call MPI_FILE_CLOSE(samplepos_filehandle, ierr)
    end if
    if (myid == 0) then
        close (99)
    end if
    call MPI_FINALIZE(ierr)
end program acceleration

subroutine start_particle(set, myid, n_proc)
    use user_variables, only: n_start!, debug
    use internal, only: n_in
    use test_var, only: n_injected, sec

    implicit none
    integer, intent(in) :: set, myid, n_proc

    n_injected = 0
    do while (n_injected < n_start)

        if (n_in == 0) then
            n_injected = n_injected + 1
            call inject !(n_injected)
            n_in = 1
            sec = 0
            !write(*,*)
            !write(*,*) 'primary',n_injected,n_start
        else
            sec = 1
            !write(*,*) 'secondary'
        end if
        call tracer(set, n_injected, n_proc)
        if (myid == 0 .and. mod(n_injected*100, n_start) == 0 .and. sec == 0) &
            write (*, *) set, n_injected*n_proc
    end do
end subroutine start_particle

subroutine tracer(set, n_injected, n_proc)
    use event_internal; use internal, only: n_in
    use random_walk

    implicit none
    integer, intent(in) :: set, n_injected, n_proc
    integer id
    integer, pointer :: A, pid

    A => event(n_in)%A
    pid => event(n_in)%pid
    if (A > 1) then
        id = 100 + A
    else
        id = abs(pid)
    end if
    select case (id)
        !case (102:108)                                       ! discard low A nuclei
    case (102:144)                                        ! discard low A nuclei
        n_in = n_in - 1
        return
    case (7, 145:159)
        call pitch_angle_random_walk(set, n_injected, n_proc)
    case default
        write (*, *) 'A,pid', A, pid
        call error('wrong particle typ in tracer', 0)
    end select
end subroutine tracer
