!=============================================================================!
!=============================================================================!
program acceleration
   use result; use user_variables, only: n_sets, outdir, num_steps_tot_str
   implicit none
   integer myid, n_proc
   integer set

! non-MPI values
   myid = 0
   n_proc = 1

   call init(myid)

   do set = 1, n_sets

      call start_particle(set, myid, n_proc)

! non-MPI values
      En_f_tot = En_f
      En_f_tot = En_f
   end do

   ! Write distance data to file
   open(20, file=trim(outdir)//'/isotropic_rw_fdist_'//trim(num_steps_tot_str), form='unformatted')
   write(20) final_distances
   close(20)
   open(20, file=trim(outdir)//'/isotropic_rw_pos_'//trim(num_steps_tot_str), form='unformatted')
   write(20) final_positions
   close(20)
   close (99)

end program acceleration
!=============================================================================!
!=============================================================================!
subroutine start_particle(set, myid, n_proc)
   use user_variables, only: n_start, debug
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
      else
         sec = 1
      end if
      call tracer(set, n_injected)
      if (myid == 0 .and. mod(n_injected*100, n_start) == 0 .and. sec == 0) &
         write (*, *) set, n_injected*n_proc

   end do

end subroutine start_particle
!=============================================================================!
!=============================================================================!
subroutine tracer(set, n_injected)
   use event_internal; use internal, only: n_in
   use user_variables, only: n_sets, n_start
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
!  case (102:108)                                        ! discard low A nuclei
   case (102:144)                                        ! discard low A nuclei
      n_in = n_in - 1
      return
   case (7, 145:159)
      print *, "Starting particle # ", n_injected + (set-1) * n_start
      call random_walk(set, n_injected)
   case default
      write (*, *) 'A,pid', A, pid
      call error('wrong particle typ in tracer', 0)
   end select

end subroutine tracer
!=============================================================================!
!=============================================================================!
subroutine random_walk(set, n_injected) ! w/wo diffusion in trapping phase
   use user_variables, only: debug, num_steps_tot
   use SNR_data, only: t_max; 
   use constants; use particle_data, only: m_p
   use event_internal; use result
   use internal
   use test_var, only: sec, accel
   implicit none
   integer k, n_step, num_steps_taken
   integer, intent(in) :: set, n_injected
   double precision r, m, f, df, dt, dE, delta, l_0, l_0_0
   double precision r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_2
   double precision ran0, R_L, t_shock, v_shock
   integer, pointer :: pid, A, Z
   double precision, pointer :: E, x(:), t, w

   num_steps_taken = 0

   pid => event(n_in)%pid
   A => event(n_in)%A
   Z => event(n_in)%Z
   E => event(n_in)%E
   x => event(n_in)%x
   t => event(n_in)%t
   w => event(n_in)%w

   r = ran0()
   m = A*m_p
   f = 0.d0

   ! Initiate at (0, 0, 0)
   x(1) = 0
   x(2) = 0
   x(3) = 0

   do 
      ! Step size
      df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)            ! interaction rate (1/yr)
      call scales_charged(m, Z, E, t, w, df, dt, dE)
      l_0 = R_L(E, t)/dble(Z)
      l_0_0 = l_0
      if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)

      ! Step direction (isotropic)
      call isotropic(phi, theta)

      ! Number of steps
      if (dt >= l_0) then
         ! one random step of size l_0
         dE = dE*l_0/dt
         dt = l_0
         n_step = 1
      else                                    
         ! n steps l0 in same direction
         print *, "more steps"
         l_0 = dt
         if (l_0_0/dt < 1.d3) then
            n_step = int(l_0_0/dt + 0.5d0)
         else                                 
            ! fast decays lead to overflow
            n_step = 1000 ! this should be enough
         end if
         if (debug > 0) write (*, *) 'E, step number', E, n_step
      end if
      if (n_step < 1) then
         write (*, *) l_0_0, l_0
         write (*, *) dt, df
         write (*, *) A, Z
         write (*, *) E
         write (*, *) l_0_0/dt, n_step
         call error('wrong step number', 0)
      end if

      ! Increment step count
      num_steps_taken = num_steps_taken + 1

      ! Perform step(s)
      do k = 1, n_step 
         d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! Old distance
       
         ! Random step
         x(1) = x(1) + l_0*cos(phi)*sin(theta)
         x(2) = x(2) + l_0*sin(phi)*sin(theta)
         x(3) = x(3) + l_0*cos(theta)

         d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! new distance
         t = t + dt
         E = E + dE

         dmax = 3.d0*l_0_0/v_2
         f = f + df*dt                      ! \int dt f(t)
         delta = exp(-f)                    ! exp(-\int dt f(t))
         
         ! Exit when max_steps reached
         if (num_steps_taken >= num_steps_tot) then
            call store(d2, set, n_injected, x(1), x(2), x(3))
            n_in = n_in - 1
            n_out = n_out + 1
            return
         end if
      end do
   end do

end subroutine random_walk
!=============================================================================!
!=============================================================================!
subroutine scales_charged(m, Z, En, t, w, df, dt, dE)
   use SNR_data, only: t_max
   use internal
   implicit none
   integer Z
   double precision m, En, t, w, df, dE, dt, tau_eff
   double precision dE_loss_syn
   double precision tausyn

   if (df > 0.d0) then
      tau_eff = 1.d0/df                               ! yr, interaction + decay
   else
      tau_eff = huge(0.d0)
   end if
   
   !!!  tausyn = tau_syn(m,En,t)/dble(Z**2)                           ! synchrotron
   tausyn = huge(0.d0)
   dt = 9.d-3*min(tau_eff, tausyn)

   dE_loss_syn = dt/tausyn*En                                             ! eV
   dE = -dE_loss_syn

   if (abs(dE)/En > 1.d-2) call error('dE too large in scale', 1)

   ! no sync. photo-production
end subroutine scales_charged
!=============================================================================!
!=============================================================================!
subroutine store(df, set, n_injected, x1, x2, x3) ! df = final distance
   use internal; use result
   use user_variables, only: n_start
   implicit none
   double precision, intent(in) :: df, x1, x2, x3
   integer, intent(in) :: set, n_injected
   integer idx
   
   idx = n_injected + (set-1)*n_start
   final_distances(idx) = df
   final_positions(1, idx) = x1
   final_positions(2, idx) = x2
   final_positions(3, idx) = x3
end subroutine store
!=============================================================================!
!=============================================================================!
