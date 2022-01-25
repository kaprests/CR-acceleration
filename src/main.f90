program acceleration
   use mpi_f08
   use result
   use user_variables, only: n_sets, filename, outdir, stepsize_exp_str

   implicit none
   integer set
   integer myid, n_proc, ierr, n_array
   type(MPI_FILE) traj_filehandle, crossang_filehandle, phase_space_filehandle
   integer traj_array_count, traj_array_bsize
   integer crossang_array_count, crossang_array_bsize
   integer phase_space_array_count, phase_space_array_bsize
   integer(kind=MPI_OFFSET_KIND) :: traj_offset, crossang_offset, phase_space_offset

   ! MPI init
   call MPI_INIT(ierr)  ! Initiate/create 
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) ! Get proc id
   call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, ierr) ! Get communicator size (# procs)

   ! Siulation init
   call init(myid, n_proc)

   ! MPI data file -- trajectories
   call init_mpi_io(MPI_COMM_WORLD, &
      traj_filehandle, &
      trim(outdir)//'/trajectories'//filename, &
      ierr &
   )
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

   ! Simulation loop -- loop through each set of particles and simulate
   do set = 1, n_sets
      call start_particle(set, myid, n_proc)

      ! Collect data from all procs
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
      call MPI_REDUCE(&
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

      ! Write trajectory data
      traj_array_bsize = sizeof(trajectories)
      traj_offset = myid * n_sets * traj_array_bsize + traj_array_bsize * (set - 1)
      traj_array_count = size(trajectories)
      call MPI_FILE_WRITE_AT(&
         traj_filehandle, &
         traj_offset, &
         trajectories, &
         traj_array_count, &
         MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, &
         ierr &
      )
      ! Write cross angle data
      crossang_array_bsize = sizeof(crossing_flight_angles)
      crossang_offset = myid * n_sets * traj_array_bsize + traj_array_bsize * (set - 1)
      crossang_array_count = size(crossing_flight_angles)
      call MPI_FILE_WRITE_AT(&
         crossang_filehandle, &
         crossang_offset, &
         crossing_flight_angles, &
         crossang_array_count, &
         MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, &
         ierr &
      )

      ! Write phase data 
      phase_space_array_bsize = sizeof(phase_space_dist)
      phase_space_offset = myid * n_sets * phase_space_array_bsize + phase_space_array_bsize * (set - 1)
      phase_space_array_count = size(phase_space_dist)
      call MPI_FILE_WRITE_AT(&
         phase_space_filehandle, &
         crossang_offset, &
         phase_space_dist, &
         phase_space_array_count, &
         MPI_DOUBLE_PRECISION, &
         MPI_STATUS_IGNORE, &
         ierr &
      )
   end do

   ! close outputs and finalize
   call MPI_FILE_CLOSE(traj_filehandle, ierr)
   call MPI_FILE_CLOSE(crossang_filehandle, ierr)
   call MPI_FILE_CLOSE(phase_space_filehandle, ierr)
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
