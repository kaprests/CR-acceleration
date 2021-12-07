program acceleration
   use mpi
   use result; use user_variables, only: n_sets, filename, outdir, stepsize_exp_str
   implicit none
   integer myid, n_proc, ierr, n_array
   integer set
   integer traj_filehandle, traj_array_count, traj_array_bsize
   integer crossang_filehandle, crossang_array_count, crossang_array_bsize
   integer(kind=MPI_OFFSET_KIND) :: traj_disp, crossang_disp

   ! non-MPI values
   !myid = 0
   !n_proc = 1

   ! MPI setup
   call MPI_INIT(ierr)  ! Initiate/create 
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) ! Get proc id
   call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr) ! Get communicator size (# procs)

   ! Simulation setup
   call init(myid)

   ! Non-MPI data files
   !open(20, &
   !   file=trim(outdir)//'/trajectories'//'_stepexp'//trim(stepsize_exp_str)//filename, &
   !   form='unformatted')
   !open(21, &
   !   file=trim(outdir)//'/cross_angles'//'_stepexp'//trim(stepsize_exp_str)//filename, &
   !   form='unformatted')

   ! MPI data file -- trajectories
   call MPI_FILE_OPEN(&
      MPI_COMM_WORLD, &
      trim(outdir)//'/trajectories'//'_stepexp'//trim(stepsize_exp_str)//filename, &
      MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
      MPI_INFO_NULL, &
      traj_filehandle, &
      ierr &
   )

   if (ierr > 0) then
      call MPI_FILE_DELETE(&
         trim(outdir)//'/trajectories'//'_stepexp'//trim(stepsize_exp_str)//filename, &
         MPI_INFO_NULL, &
         ierr &
      )

      call MPI_FILE_OPEN(&
         MPI_COMM_WORLD, &
         trim(outdir)//'/trajectories'//'_stepexp'//trim(stepsize_exp_str)//filename, &
         MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
         MPI_INFO_NULL, &
         traj_filehandle, &
         ierr &
      )
   end if

   ! MPI data file -- crossing flight angles
   call MPI_FILE_OPEN(&
      MPI_COMM_WORLD, &
      trim(outdir)//'/cross_angles'//'_stepexp'//trim(stepsize_exp_str)//filename, &
      MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
      MPI_INFO_NULL, &
      crossang_filehandle, &
      ierr &
   )

   if (ierr > 0) then
      call MPI_FILE_DELETE(&
         trim(outdir)//'/cross_angles'//'_stepexp'//trim(stepsize_exp_str)//filename, &
         MPI_INFO_NULL, &
         ierr &
      )

      call MPI_FILE_OPEN(&
         MPI_COMM_WORLD, &
         trim(outdir)//'/cross_angles'//'_stepexp'//trim(stepsize_exp_str)//filename, &
         MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
         MPI_INFO_NULL, &
         crossang_filehandle, &
         ierr &
      )
   end if

   do set = 1, n_sets

      call start_particle(set, myid, n_proc)

      ! non-MPI values
      !En_f_tot = En_f

      n_array = (2*pid_max+1)*n_enbin

      call MPI_REDUCE(&
         En_f, &
         En_f_tot, &
         n_array, &
         MPI_DOUBLE_PRECISION, &
         MPI_SUM, &
         0, &
         MPI_COMM_WORLD, &
         ierr &
      )

      if (myid == 0) call output(set, n_proc)

      ! Non-MPI write
      !write(20) trajectories
      !write(21) crossing_flight_angles

      ! Write trajectory data
      traj_array_bsize = sizeof(trajectories)
      traj_disp = traj_array_bsize * set * myid
      traj_array_count = size(trajectories)

      call MPI_FILE_SET_VIEW(&
         traj_filehandle, &
         traj_disp, &
         MPI_INTEGER, &
         MPI_INTEGER, &
         'native', &
         MPI_INFO_NULL, &
         ierr &
      )
      call MPI_FILE_WRITE(&
         traj_filehandle, &
         trajectories, &
         traj_array_count, &
         MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, &
         ierr &
      )

      ! Write cross angle data
      crossang_array_bsize = sizeof(crossing_flight_angles)
      crossang_disp = traj_array_bsize * myid
      crossang_array_count = size(crossing_flight_angles)

      call MPI_FILE_SET_VIEW(&
         crossang_filehandle, &
         crossang_disp, &
         MPI_INTEGER, &
         MPI_INTEGER, &
         'native', &
         MPI_INFO_NULL, &
         ierr &
      )
      call MPI_FILE_WRITE(&
         crossang_filehandle, &
         crossing_flight_angles, &
         crossang_array_count, &
         MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, &
         ierr &
      )

   end do

   !close(20)
   !close(21)
   !call output_finish

   call MPI_FILE_CLOSE(traj_filehandle, ierr)
   call MPI_FILE_CLOSE(crossang_filehandle, ierr)

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
   use acceleration
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
      if (isotropic) then
         call isotropic_random_walk(set, n_injected)
      else
         call pitch_angle_accel(set, n_injected)
      end if
   case default
      write (*, *) 'A,pid', A, pid
      call error('wrong particle typ in tracer', 0)
   end select

end subroutine tracer
