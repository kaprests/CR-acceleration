!=============================================================================!
!=============================================================================!
program acceleration
   use result; use user_variables, only: n_sets, outdir, theta_max_str, t_max_str
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
   open(20, file=trim(outdir)//'/pitch_angle_rw_fdist_tmax'//trim(t_max_str)//'_theta'//theta_max_str, form='unformatted')
   write(20) final_distances
   close(20)
   open(20, file=trim(outdir)//'/pitch_angle_rw_pos_tmax'//trim(t_max_str)//'_theta'//theta_max_str, form='unformatted')
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
   use user_variables, only: debug, t_max
   use constants; use particle_data, only: m_p
   use event_internal; use result
   use internal
   use test_var, only: sec, accel

   implicit none

   integer n_step, num_steps_taken
   integer, intent(in) :: set, n_injected
   double precision r, m, f, df, dt, dE, delta, l_0, l_0_0
   double precision r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_2
   double precision ran0, R_L, t_shock, v_shock, get_v_2
   integer, pointer :: pid, A, Z
   double precision, pointer :: E, x(:), t, w
   double precision :: g(3), p(3), R_euler(3,3), theta_e, phi_e
   integer i, j, k
   double precision :: gamma_v, cos_theta, theta_max_computed

   num_steps_taken = 0

   pid => event(n_in)%pid
   A => event(n_in)%A
   Z => event(n_in)%Z
   E => event(n_in)%E
   x => event(n_in)%x
   t => event(n_in)%t
   w => event(n_in)%w

   d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)    
   if (sec == 0 .and. abs(d1/t_shock(t) - 1.01d0) .gt. 1.d-6) then    
      call error('wrong initial condition, shock', 0)    
   end if

   r = ran0()
   m = A*m_p
   f = 0.d0

   do 
      ! Step size
      df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)            ! interaction rate (1/yr)
      call scales_charged(m, Z, E, t, w, df, dt, dE)
      l_0 = R_L(E, t)/dble(Z)

      ! Modify stepsize for pitch angle scattering
      call max_scattering_angle(theta_max_computed, v_shock(t), E)
      l_0 = l_0 * (theta_max_computed/pi)**2 ! Guess work
      l_0_0 = l_0 ! Why l_0_0
      if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)

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
         call scattering_angle(theta_e, phi_e, theta_max_computed)

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

      if (num_steps_taken == 1000) then
         print *, "exceeded 1000 steps"
      end if

      if (num_steps_taken == 10000) then
         print *, "exceeded 10000 steps"
      end if


      ! Perform step(s)
      do k = 1, n_step 
         d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! radial dist before step
         r_sh1 = t_shock(t)                     ! shock radius before step
       
         ! Random step
         if (d1 < r_sh1) then
            ! Particle in dowstream -> advection
            call radially_outward(phi_v, theta_v, x(1), x(2), x(3))
            v_2 = get_v_2(v_shock(t)) ! approaching velocity of US seen by DS
            if (v_2 <= 0.d0) call error('v_2<=0', 0)
            
            ! Lorentz transformed step
            gamma_v = 1.d0/sqrt(1.d0 - v_2**2) ! ([v_2]=1 -> v_2 = beta = v/c)
            x(1) = x(1) + v_2*cos(phi_v)*sin(theta_v)*dt + l_0*sin(theta)*cos(phi)/gamma_v    
            x(2) = x(2) + v_2*sin(phi_v)*sin(theta_v)*dt + l_0*sin(theta)*sin(phi)/gamma_v    
            x(3) = x(3) + v_2*cos(theta_v)*dt + l_0*cos(theta)/gamma_v
         else 
            ! Particle in upstream -> 'normal' step
            x(1) = x(1) + l_0*cos(phi)*sin(theta)
            x(2) = x(2) + l_0*sin(phi)*sin(theta)
            x(3) = x(3) + l_0*cos(theta)
         end if

         t = t + dt
         E = E + dE
         d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2) ! radial dist after step
         r_sh2 = t_shock(t)                     ! shock radius after step

         ! Check for shock crossing
         if (d2 < r_sh2 .and. r_sh1 < d1) then
            ! Crossed shock (US -> DS)
            call radially_outward(phi_v, theta_v, x(1), x(2), x(3))
            v_2 = get_v_2(v_shock(t)) ! US sees DS approach at velocity v_2
            if (v_2 <= 0.d0) call error('v_2<=0', 0)    

            cos_theta = &    
               cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &    
               sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &    
               cos(theta)*cos(theta_v)    

            gamma_v = 1.d0/sqrt(1.d0 - v_2**2)                 ! v_2 dimless (v_2 = beta = v/c)    
            !E_old = E                                          ! Old energy    
            E = gamma_v*E*(1 - v_2*cos_theta)                  ! New energy
         else if (d2 > r_sh2 .and. r_sh1 > d1) then
            ! Crossed shock (DS -> US)
            call radially_outward(phi_v, theta_v, x(1), x(2), x(3)) ! direction of v_2 and shock
            v_2 = get_v_2(v_shock(t)) ! DS sees US approach at same velocity v_2
            if (v_2 <= 0.d0) call error('v_2<=0', 0)

            cos_theta = &
               cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
               sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
               cos(theta)*cos(theta_v)

            gamma_v = 1.d0/sqrt(1.d0 - v_2**2)
            !E_old = E
            E = gamma_v*E*(1 + v_2*cos_theta)
         end if

         dmax = 3.d0*l_0_0/v_2
         f = f + df*dt                      ! \int dt f(t)
         delta = exp(-f)                    ! exp(-\int dt f(t))
         
         ! Exit when t_max exceeded
         if (t > t_max .or. d2 < r_sh2 - dmax .or. r > delta) then
            if (t > t_max .or. d2 < r_sh2 - dmax) then
               !print *, "t max: ", t_max
               !print *, "steps taken: ", num_steps_taken
               !call store(d2, set, n_injected, x(1), x(2), x(3))
               call store(pid, E, w)
               n_in = n_in - 1
               n_out = n_out + 1
               return
            end if
            if (r > delta) then
               write (*, *) 'should never happen...'
               stop
            end if
         end if
      end do
   end do

end subroutine random_walk
!=============================================================================!
!=============================================================================!
subroutine scales_charged(m, Z, En, t, w, df, dt, dE)
   use user_variables, only: t_max
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
subroutine store(pid, En, w) !, num_crossings, rel_energy_gain_sum)
   use internal; use result, only : n_enbin, En_f, NE_esc, rel_energy_gain_total_sum
   implicit none
   integer pid, i
   double precision En, w, l
!   integer, intent(in) :: num_crossings
!   double precision, intent(in) :: rel_energy_gain_sum
   double precision  :: rel_energy_gain_avg

   l = log10(En)                                                ! energy bin
   i = int((l - d_f)/dn) ! dn = 0.1d0

   if (i <= 0) then
      call error('stor_esc, E<=Emin', 11)
      i = 1
   end if
   if (i > n_enbin) then
      i = n_enbin
   end if

   En_f(pid, i) = En_f(pid, i) + w*En
   NE_esc(i) = NE_esc(i) + 1
!   rel_energy_gain_avg = rel_energy_gain_sum / num_crossings
!   rel_energy_gain_total_sum = rel_energy_gain_total_sum + rel_energy_gain_sum
!  write(*,*) 'store: ',pid,i
end subroutine store
!=============================================================================!
!=============================================================================!
