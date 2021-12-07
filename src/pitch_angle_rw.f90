program pitch_angle
   use result; use user_variables, only: n_sets, outdir, theta_max_str, t_max_str, stepsize_exp, stepsize_exp_str
   use internal, only: E_inj
   implicit none
   integer myid, n_proc, set

   E_inj = 1.0d11

   ! non-MPI values
   myid = 0
   n_proc = 1

   call init(myid)
   
   open(20, &
      file=trim(outdir)//'/pitch_angle_rw_pos_tmax'//trim(t_max_str)//'_theta'//&
      trim(theta_max_str)//'_stepexp'//trim(stepsize_exp_str),form='unformatted')
   open(21, &
      file=trim(outdir)//'/pitch_angle_rw_trajectories_tmax'//&
      trim(t_max_str)//'_theta'//trim(theta_max_str)//'_stepexp'//trim(stepsize_exp_str), &
      form='unformatted')
   open(22, &
      file=trim(outdir)//'/pitch_angle_rw_samplepos_tmax'//&
      trim(t_max_str)//'_theta'//trim(theta_max_str)//'_stepexp'//trim(stepsize_exp_str), &
      form='unformatted')

   do set = 1, n_sets

      call start_particle(set, myid, n_proc)

      ! non-MPI values
      En_f_tot = En_f
      En_f_tot = En_f

      write(20) final_positions
      write(21) trajectories
      write(22) sample_positions
   end do

   close(20)
   close(21)
   close(22)

   close (99)

end program pitch_angle


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
      else
         sec = 1
      end if
      call tracer(set, n_injected)
      if (myid == 0 .and. mod(n_injected*100, n_start) == 0 .and. sec == 0) &
         write (*, *) set, n_injected*n_proc

   end do

end subroutine start_particle


subroutine tracer(set, n_injected)
   use event_internal; use internal, only: n_in
   use user_variables, only: n_start
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


subroutine random_walk(set, n_injected) ! w/wo diffusion in trapping phase
   use user_variables, only: debug, t_max, theta_max, num_steps_log, stepsize_exp
   use constants; use particle_data, only: m_p
   use event_internal; use result
   use internal
   implicit none
   integer n_step, num_steps_taken
   integer, intent(in) :: set, n_injected
   double precision r, m, f, df, dt, dE, delta, l_0, l_0_0
   double precision phi, theta, d1, d2, dmax, v_2
   double precision ran0, R_L
   integer, pointer :: pid, A, Z
   double precision, pointer :: E, x(:), t, w
   double precision :: g(3), p(3), R_euler(3,3), theta_e, phi_e
   integer i, j, k
   double precision :: t0
   integer :: num_steps_total, idx, sample_count, sample_int, num_samples

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

   ! Step size
   df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)            ! interaction rate (1/yr)
   call scales_charged(m, Z, E, t, w, df, dt, dE)
   l_0 = R_L(E, t)/dble(Z)

   ! Modify stepsize for pitch angle scattering
   l_0 = l_0 * (theta_max/pi)**stepsize_exp ! Guess work

   l_0_0 = l_0 ! Why l_0_0
   if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)

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

   t0 = t
   num_steps_total = abs((t0-t_max))/dt+1
   sample_int = floor(real(num_steps_total/num_steps_log)) + 1 ! sample interval
   num_samples = floor(real(num_steps_total/sample_int))
   sample_count = 0
   num_steps_taken = 0
   if (.not. allocated(sample_positions)) allocate(sample_positions(4, n_start, num_samples))

   do 
      ! log position
      if (num_steps_taken+1 <= size(trajectories, 3)) then
         trajectories(1, n_injected, num_steps_taken+1) = x(1)
         trajectories(2, n_injected, num_steps_taken+1) = x(2)
         trajectories(3, n_injected, num_steps_taken+1) = x(3)
         trajectories(4, n_injected, num_steps_taken+1) = t
      end if

      ! sample positions
      if (modulo(num_steps_taken+1,sample_int)==0) then
         if (sample_count > num_samples) then
            call error('Sample count to big!', 1)
         end if
         if (sample_count == num_samples) then
            print *, "Sample filled!"
         end if
         sample_count = sample_count + 1
         sample_positions(1, n_injected, sample_count) = x(1)
         sample_positions(2, n_injected, sample_count) = x(2)
         sample_positions(3, n_injected, sample_count) = x(3)
         sample_positions(4, n_injected, sample_count) = t
      end if

      ! Step direction 
      if (num_steps_taken == 0) then
         ! First step isotropic
         call isotropic(phi, theta)
      else
         ! Then small angle
         ! Initial momentum vector
         g(1) = cos(phi)*sin(theta)    
         g(2) = sin(phi)*sin(theta)    
         g(3) = cos(theta)

         ! scattering angles in rotated frame, where p=(1,0,0)
         call scattering_angle_dev(theta_e, phi_e)

         ! Scattered momentum vector in rotated frame
         p(1) = cos(phi_e)*sin(theta_e)
         p(2) = sin(phi_e)*sin(theta_e)    
         p(3) = cos(theta_e)

         ! Rotate back to lab frame coordinates
         call euler_RyRz(-theta, -phi, R_euler)
         g = 0.d0                                      
         do i=1,3 ! for column i in columns
            do j=1,3 ! for row j in rows       
               g(i) = g(i) + R_euler(i,j)*p(j)
            end do                                
         end do

         ! New theta and phi from cartesian coordinates
         theta = atan2( sqrt(g(1)**2+g(2)**2), g(3) )              
         phi   = atan( g(2)/g(1) )                       
         if (g(1)<0.d0.and.g(2)>0) phi = phi+pi
         if (g(1)<0.d0.and.g(2)<0) phi = phi+pi    
         if (g(1)>0.d0.and.g(2)<0) phi = phi+two_pi
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
         
         ! Exit when t_max exceeded
         if (t > t_max) then
            idx = n_injected + (set-1)*n_start
            if (idx == n_start*n_sets) then
               print *, "t max: ", t_max
               print *, "steps taken: ", num_steps_taken
               print *, "num_steps_total: ", num_steps_total
               print *, "num_steps_log: ", num_steps_log
               print *, "sample_count: ", sample_count
               print *, "num_samples: ", num_samples
            end if
            call store(n_injected, x(1), x(2), x(3))
            n_in = n_in - 1
            n_out = n_out + 1
            return
         end if
      end do
   end do

end subroutine random_walk


subroutine scales_charged(m, Z, En, t, w, df, dt, dE)
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


subroutine store(n_injected, x1, x2, x3)
   use result, only: final_positions
   implicit none
   double precision, intent(in) :: x1, x2, x3
   integer, intent(in) :: n_injected

   final_positions(1, n_injected) = x1
   final_positions(2, n_injected) = x2
   final_positions(3, n_injected) = x3

end subroutine store