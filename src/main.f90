! File: main.f90
! Cosmic ray acceleration

program acceleration
    use mpi
    use result; use user_variables, only: n_sets

    implicit none
    integer :: myid, n_proc, ierr, n_array
    integer :: set

    ! non-MPI values
    !myid = 0
    !n_proc = 1

    ! MPI init
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr)

    ! Initialize simulation
    call init(myid, n_proc)

    do set = 1, n_sets
        call start_particle(set, myid, n_proc)

        ! non-MPI values
        !En_f_tot = En_f

        n_array = (2*pid_max + 1)*n_enbin
        call MPI_REDUCE( &
            En_f, En_f_tot, n_array, MPI_DOUBLE_PRECISION, MPI_SUM, 0, &
            MPI_COMM_WORLD, ierr &
            )

        if (myid == 0) call output(set, n_proc)
    end do
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
        call tracer(set, n_injected)
        if (myid == 0 .and. mod(n_injected*100, n_start) == 0 .and. sec == 0) &
            write (*, *) set, n_injected*n_proc
    end do
end subroutine start_particle

subroutine tracer(set, n_injected)
    use event_internal; use internal, only: n_in
    use random_walk

    implicit none
    integer, intent(in) :: set, n_injected
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
        call pitch_angle_random_walk(set, n_injected)
    case default
        write (*, *) 'A,pid', A, pid
        call error('wrong particle typ in tracer', 0)
    end select
end subroutine tracer
