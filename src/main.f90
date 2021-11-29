program acceleration
   use result; use user_variables, only: n_sets, filename, outdir, stepsize_exp_str
   implicit none
   integer myid, n_proc !,ierr,n_array
   integer set

   ! non-MPI values
   myid = 0
   n_proc = 1

   call init(myid)

   open(20, &
      file=trim(outdir)//'/trajectories'//'_stepexp'//trim(stepsize_exp_str)//filename, &
      form='unformatted')

   do set = 1, n_sets

      call start_particle(set, myid, n_proc)

      ! non-MPI values
      En_f_tot = En_f
      En_f_tot = En_f

      if (myid == 0) call output(set, n_proc)

      write(20) trajectories
   end do

   close(20)
   !call output_finish

   close (99)

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
