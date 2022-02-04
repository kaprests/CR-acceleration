module acceleration
   implicit none
   private

   public isotropic_random_walk, pitch_angle_accel
contains

   subroutine pitch_angle_accel(set, n_injected)
      ! Pitch angle scattering shock acceleration
      use user_variables, only: debug, t_max, num_steps_log, stepsize_exp, inj_model
      use constants
      use particle_data, only: m_p
      use event_internal
      use result
      use internal
      use test_var, only: sec, accel

      implicit none
      ! Dummy variables
      integer, intent(in) :: set, n_injected
      ! Loop variables: incr, end
      integer :: k, n_step
      ! Functions
      double precision :: ran0, R_L, t_shock, v_shock, get_v_2, cubic_spline_stepsize, v_particle
      double precision :: analytical_stepsize
      ! Pointers
      integer, pointer :: pid, A, Z
      double precision, pointer :: E, x(:), t, w
      ! Variables
      double precision :: r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_2, v_p
      double precision :: r, m, f, df, dt, dE, delta, l_0, l_0_0
      double precision :: v_x, v_y, v_z
      double precision :: v_px, v_py, v_pz ! particle velocity components
      double precision :: gamma_v, gamma_x, gamma_y, gamma_z, cos_theta
      double precision :: gamma_pv, gamma_px, gamma_py, gamma_pz
      double precision :: l_x, l_y, l_z
      double precision :: d0, r_sh0
      ! ############
      integer :: num_crossings, num_steps_taken
      double precision :: rel_energy_gain, E_old, rel_energy_gain_sum
      double precision :: theta_max, theta_max0
      double precision :: g(3), p(3), R_euler(3, 3), theta_e, phi_e
      integer :: i, j, idx

      pid => event(n_in)%pid
      A => event(n_in)%A
      Z => event(n_in)%Z
      E => event(n_in)%E
      x => event(n_in)%x
      t => event(n_in)%t
      w => event(n_in)%w

      ! Initial radial position of particle
      d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
      if (sec == 0 .and. abs(d1/t_shock(t) - 1.01d0) .gt. 1.d-6 .and. inj_model == 0) then
         call error('wrong initial condition, shock', 0)
      else if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6 .and. inj_model == 1) then
         call error('wrong initial condition, shock', 0)
      else if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6 .and. inj_model == 2) then
         call error('wrong initial condition, shock', 0)
      end if

      r = ran0()
      m = A*m_p
      f = 0.d0

      rel_energy_gain_sum = 0
      num_crossings = 0
      num_steps_taken = 0

      do
         ! Max scattering angle 
         call max_scattering_angle(theta_max, v_shock(t), E)
         if (num_steps_taken == 0) then
            print *, "theta_max (initial): ", theta_max
            print *, "E_inj_exp: ", log10(E_inj)
            ! Log first theta max 
            theta_max0 = theta_max
         end if

         ! Stepsize
         df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)   ! interaction rate (1/yr)
         call scales_charged(m, Z, E, t, w, df, dt, dE)
         !l_0 = cubic_spline_stepsize(theta_max)/dble(Z)
         l_0 = analytical_stepsize(E, t, theta_max)/dble(Z)
         l_0_0 = l_0
         ! Old adjustment of stepsize
         !l_0 = l_0*(theta_max/pi)**stepsize_exp
         !l_0_0 = l_0
         if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)

         ! Step direction
         r_sh0 = t_shock(t)                              ! shock position (not needed here?)
         d0 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! particle radial position
         if (num_steps_taken == 0) then
            ! First step isotropic
            call isotropic(phi, theta)

            ! First step in radial direction
            !call radially_outward(phi, theta, x(1), x(2), x(3))
         else
            ! Following steps small angle
            ! g: initial momentum vector
            g(1) = cos(phi)*sin(theta)
            g(2) = sin(phi)*sin(theta)
            g(3) = cos(theta)

            ! Random scattering angles within scattering cone in rotated frame
            call scattering_angle(theta_e, phi_e, theta_max)

            ! Scattered momentum vector in rotated frame
            p(1) = cos(phi_e)*sin(theta_e)
            p(2) = sin(phi_e)*sin(theta_e)
            p(3) = cos(theta_e)

            ! Rotate p back to lab frame and update momentum vector
            call euler_RyRz(-theta, -phi, R_euler)
            g = 0.d0
            do i = 1, 3, 1
               do j = 1, 3, 1
                  g(i) = g(i) + R_euler(i, j)*p(j)
               end do
            end do

            ! New angles
            theta = atan2(sqrt(g(1)**2 + g(2)**2), g(3))
            phi = atan(g(2)/g(1))
            if (g(1) < 0.d0 .and. g(2) > 0) phi = phi + pi
            if (g(1) < 0.d0 .and. g(2) < 0) phi = phi + pi
            if (g(1) > 0.d0 .and. g(2) < 0) phi = phi + two_pi
         end if

         ! Number of steps
         if (dt >= l_0) then                     ! one random step of size l_0
            dE = dE*l_0/dt
            dt = l_0
            n_step = 1
         else                                    ! n steps l0 in same direction
            l_0 = dt
            if (l_0_0/dt < 1.d3) then
               n_step = int(l_0_0/dt + 0.5d0)
            else                                 ! fast decays lead to overflow
               n_step = 1000                     ! this should be enough
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

         ! Increment number of steps taken
         num_steps_taken = num_steps_taken + 1

         ! Log position
         if (num_steps_taken + 1 <= size(trajectories, 3)) then
            trajectories(1, n_injected, num_steps_taken + 1) = t
            trajectories(2, n_injected, num_steps_taken + 1) = x(1)
            trajectories(3, n_injected, num_steps_taken + 1) = x(2)
            trajectories(4, n_injected, num_steps_taken + 1) = x(3)
         end if

         if (num_steps_taken + 1 <= size(phase_space_dist, 3)) then
            ! Store in shock rest frame -- eventually transform afterwards
            v_p = v_particle(E, t)
            v_px = v_p*cos(phi)*sin(theta)
            v_py = v_p*sin(phi)*sin(theta)
            v_pz = v_p*cos(theta)
            gamma_px = 1.d0/sqrt(1.d0 - v_px**2)                 ! v_2 dimless (v_2 = beta = v/c)
            gamma_py = 1.d0/sqrt(1.d0 - v_py**2)                 ! v_2 dimless (v_2 = beta = v/c)
            gamma_pz = 1.d0/sqrt(1.d0 - v_pz**2)                 ! v_2 dimless (v_2 = beta = v/c)

            v_2 = v_shock(t)
            v_x = v_p*cos(phi_v)*sin(theta_v)
            v_y = v_p*sin(phi_v)*sin(theta_v)
            v_z = v_p*cos(theta_v)
            gamma_x = 1.d0/sqrt(1.d0 - v_x**2)                 ! v_2 dimless (v_2 = beta = v/c)
            gamma_y = 1.d0/sqrt(1.d0 - v_y**2)                 ! v_2 dimless (v_2 = beta = v/c)
            gamma_z = 1.d0/sqrt(1.d0 - v_z**2)                 ! v_2 dimless (v_2 = beta = v/c)
            gamma_v = 1.d0/sqrt(1.d0 - v_2**2)

            phase_space_dist(1, n_injected, num_steps_taken + 1) = &
               gamma_v*(t - v_2*sqrt(x(1)**2+x(2)**2+x(3)**2))
            phase_space_dist(2, n_injected, num_steps_taken + 1) = gamma_x*(x(1)-v_x*t)
            phase_space_dist(3, n_injected, num_steps_taken + 1) = gamma_y*(x(2)-v_y*t)
            phase_space_dist(4, n_injected, num_steps_taken + 1) = gamma_z*(x(3)-v_z*t)

            phase_space_dist(5, n_injected, num_steps_taken + 1) = &
               gamma_x*(gamma_px * m_p * v_px - E*v_x)
            phase_space_dist(6, n_injected, num_steps_taken + 1) = &
               gamma_y*(gamma_py * m_p * v_py - E*v_y)
            phase_space_dist(7, n_injected, num_steps_taken + 1) = &
               gamma_z*(gamma_pz * m_p * v_pz - E*v_z)
         end if
         !print *, "l_0: ", l_0
         !print *, "v: ", l_0/dt
         !print *, "D: ", l_0*(l_0/dt)/6
         !print *, "D': ", l_0*(l_0/dt)/3

         ! Perform step(s)
         do k = 1, n_step
            ! distances before step
            r_sh1 = t_shock(t)                              ! old shock position
            d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! old distance/particle radial position

            ! Perform random step
            if (d1 < r_sh1) then ! Particle in downstream
               ! find direction of v_2 from position x (radially outward):
               call radially_outward(phi_v, theta_v, x(1), x(2), x(3))

               ! v_2: velocity of downstream as seen in US
               v_2 = get_v_2(v_shock(t))
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               ! Lorentz transformed random step (advection)
               v_x = v_2*cos(phi_v)*sin(theta_v)
               v_y = v_2*sin(phi_v)*sin(theta_v)
               v_z = v_2*cos(theta_v)

               gamma_x = 1.d0/sqrt(1.d0 - v_x**2)                 ! v_2 dimless (v_2 = beta = v/c)
               gamma_y = 1.d0/sqrt(1.d0 - v_y**2)                 ! v_2 dimless (v_2 = beta = v/c)
               gamma_z = 1.d0/sqrt(1.d0 - v_z**2)                 ! v_2 dimless (v_2 = beta = v/c)

               x(1) = x(1) + v_x*dt + l_0*sin(theta)*cos(phi)/gamma_x
               x(2) = x(2) + v_y*dt + l_0*sin(theta)*sin(phi)/gamma_y
               x(3) = x(3) + v_z*dt + l_0*cos(theta)/gamma_z
            else ! Particle in upstream
               ! random step (isotropic in current medium rest frame)
               x(1) = x(1) + l_0*cos(phi)*sin(theta)
               x(2) = x(2) + l_0*sin(phi)*sin(theta)
               x(3) = x(3) + l_0*cos(theta)
            end if

            ! distances after step
            d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)              ! new particle distance
            t = t + dt
            E = E + dE
            r_sh2 = t_shock(t)                                  ! new shock dist

            ! Check for shock crossing
            if (d2 < r_sh2 .and. r_sh1 < d1) then ! Crossed shock (US->DS)
               call radially_outward(phi_v, theta_v, x(1), x(2), x(3)) ! angle of v_2 at pos x
               v_2 = get_v_2(v_shock(t)) ! US sees DS approach at velocity v_2
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               cos_theta = &
                  cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                  sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                  cos(theta)*cos(theta_v)

               gamma_v = 1.d0/sqrt(1.d0 - v_2**2)                 ! v_2 dimless (v_2 = beta = v/c)
               E_old = E                                          ! Old energy
               E = gamma_v*E*(1 - v_2*cos_theta)                  ! New energy
               rel_energy_gain = (E - E_old)/E_old                ! rel energy gain
               rel_energy_gain_sum = rel_energy_gain_sum + rel_energy_gain
               accel = 1
               num_crossings = num_crossings + 1

               ! log angles at crossing
               if (num_crossings <= size(crossing_flight_angles, 2)) then
                  crossing_flight_angles(n_injected, num_crossings) = acos(cos_theta)
               end if

               ! Lorentz tranform three momentum 
               v_x = v_2*cos(phi_v)*sin(theta_v)      ! Shock velocity in x direction
               v_y = v_2*sin(phi_v)*sin(theta_v)      ! Shock velocity in y direction
               v_z = v_2*cos(theta_v)                 ! Shock velocity in z direction

               gamma_x = 1.d0/sqrt(1.d0 - v_x**2)     ! v_2 dimless (v_2 = beta = v/c)
               gamma_y = 1.d0/sqrt(1.d0 - v_y**2)     ! v_2 dimless (v_2 = beta = v/c)
               gamma_z = 1.d0/sqrt(1.d0 - v_z**2)     ! v_2 dimless (v_2 = beta = v/c)

               !l_x = -v_x*dt + l_0*sin(theta)*cos(phi)/gamma_x
               !l_y = -v_y*dt + l_0*sin(theta)*sin(phi)/gamma_y
               !l_z = -v_z*dt + l_0*cos(theta)/gamma_z
               !call radially_outward(phi, theta, l_x, l_y, l_z)
            else if (d2 > r_sh2 .and. r_sh1 > d1) then ! Cossed shock (DS -> US)
               call radially_outward(phi_v, theta_v, x(1), x(2), x(3)) ! direction of v_2 and shock
               v_2 = get_v_2(v_shock(t)) ! DS sees US approach at same velocity v_2
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               cos_theta = &
                  cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                  sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                  cos(theta)*cos(theta_v)

               gamma_v = 1.d0/sqrt(1.d0 - v_2**2)
               E_old = E
               E = gamma_v*E*(1 + v_2*cos_theta)
               rel_energy_gain = (E - E_old)/E_old
               rel_energy_gain_sum = rel_energy_gain_sum + rel_energy_gain
               accel = 1
               num_crossings = num_crossings + 1

               if (num_crossings <= size(crossing_flight_angles, 2)) then
                  crossing_flight_angles(n_injected, num_crossings) = acos(cos_theta)
               end if

               ! Lorentz tranform three momentum 
               v_x = v_2*cos(phi_v)*sin(theta_v)
               v_y = v_2*sin(phi_v)*sin(theta_v)
               v_z = v_2*cos(theta_v)

               gamma_x = 1.d0/sqrt(1.d0 - v_x**2)                 ! v_2 dimless (v_2 = beta = v/c)
               gamma_y = 1.d0/sqrt(1.d0 - v_y**2)                 ! v_2 dimless (v_2 = beta = v/c)
               gamma_z = 1.d0/sqrt(1.d0 - v_z**2)                 ! v_2 dimless (v_2 = beta = v/c)

               !l_x = v_x*dt + l_0*sin(theta)*cos(phi)/gamma_x
               !l_y = v_y*dt + l_0*sin(theta)*sin(phi)/gamma_y
               !l_z = v_z*dt + l_0*cos(theta)/gamma_z
               !call radially_outward(phi, theta, l_x, l_y, l_z)
            end if

            v_2 = get_v_2(v_shock(t))
            dmax = 3.d0*l_0_0/v_2
            f = f + df*dt                      ! \int dt f(t)
            delta = exp(-f)                    ! exp(-\int dt f(t))

            ! exit, if a) too late, b) too far down-stream, or c) scattering:
            if (t > t_max .or. d2 < r_sh2 - dmax .or. r > delta) then
               if (t > t_max .or. d2 < r_sh2 - dmax) then  ! we're tired or trapped behind
                  ! write(*,*) 'tired',n_in,n_out
                  r_sh0 = t_shock(t)
                  d0 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
                  if (t > t_max .and. d0 < r_sh0) then
                     print *, "Time exit - particle in downstream - num_cross: ", num_crossings
                  else if (t > t_max .and. d0 > r_sh0) then
                     print *, "Time exit - particle in upstream - num_cross: ", num_crossings
                  end if
                  print *, "!!!!!!!!!!!!!!!!!"
                  print *, "w: ", w
                  print *, "!!!!!!!!!!!!!!!!!"
                  call store(pid, E, w, num_crossings, rel_energy_gain_sum)
                  !call store_raw(E, set, n_injected, num_crossings)
                  print *, "#############################"
                  print *, "Num shock crossings: ", num_crossings
                  print *, "Num steps taken: ", num_steps_taken
                  print *, "Initial theta max: ", theta_max0
                  print *, "Final theta max: ", theta_max
                  n_in = n_in - 1
                  n_out = n_out + 1
                  return
               end if
               if (r > delta) then                              ! decay or scattering
                  write (*, *) 'should never happen...'
                  stop
               end if
            end if

         end do
      end do
   end subroutine pitch_angle_accel

   subroutine isotropic_random_walk(set, n_injected)
      ! Isotropic random walk shock acceleration
      use user_variables, only: debug, t_max, inj_model
      use constants; 
      use particle_data, only: m_p
      use event_internal
      use result
      use internal
      use test_var, only: sec, accel

      implicit none

      integer, intent(in) :: set, n_injected
      integer k, n_step
      double precision r, m, f, df, dt, dE, delta, l_0, l_0_0
      double precision r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_2
      double precision ran0, R_L, t_shock, v_shock, get_v_2
      integer, pointer :: pid, A, Z
      double precision, pointer :: E, x(:), t, w
      double precision :: v_x, v_y, v_z ! Shock velocity components
      double precision :: gamma_v, gamma_x, gamma_y, gamma_z, cos_theta
      double precision :: d0, r_sh0

      ! ############
      integer :: num_crossings, num_steps_taken
      double precision :: rel_energy_gain, E_old, rel_energy_gain_sum

      pid => event(n_in)%pid
      A => event(n_in)%A
      Z => event(n_in)%Z
      E => event(n_in)%E
      x => event(n_in)%x
      t => event(n_in)%t
      w => event(n_in)%w

      d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
      if (sec == 0 .and. abs(d1/t_shock(t) - 1.01d0) .gt. 1.d-6 .and. inj_model == 0) then
         call error('wrong initial condition, shock', 0)
      else if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6 .and. inj_model == 1) then
         call error('wrong initial condition, shock', 0)
      else if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6 .and. inj_model == 2) then
         call error('wrong initial condition, shock', 0)
      end if

      r = ran0()
      m = A*m_p
      f = 0.d0

      rel_energy_gain_sum = 0
      num_crossings = 0
      num_steps_taken = 0

      do
         ! Find step size
         df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)   ! interaction rate (1/yr)
         call scales_charged(m, Z, E, t, w, df, dt, dE)
         l_0 = R_L(E, t)/dble(Z)
         l_0_0 = l_0
         if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)

         ! Find step direction
         r_sh0 = t_shock(t)                              ! shock position
         d0 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! particle radial position
         ! Non relativistic -> isotropic random walk
         call isotropic(phi, theta)

         ! Number of steps
         if (dt >= l_0) then                     ! one random step of size l_0
            dE = dE*l_0/dt
            dt = l_0
            n_step = 1
         else                                    ! n steps l0 in same direction
            l_0 = dt
            if (l_0_0/dt < 1.d3) then
               n_step = int(l_0_0/dt + 0.5d0)
            else                                 ! fast decays lead to overflow
               n_step = 1000                     ! this should be enough
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

         ! Increment steps taken
         num_steps_taken = num_steps_taken + 1

         ! Log position
         if (num_steps_taken + 1 <= size(trajectories, 3)) then
            trajectories(1, n_injected, num_steps_taken + 1) = x(1)
            trajectories(2, n_injected, num_steps_taken + 1) = x(2)
            trajectories(3, n_injected, num_steps_taken + 1) = x(3)
            trajectories(4, n_injected, num_steps_taken + 1) = t
         end if

         ! perform random step(s)
         do k = 1, n_step
            ! distances before step
            r_sh1 = t_shock(t)                              ! old shock position
            d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! old distance/particle radial position

            ! Perform random step
            if (d1 < r_sh1) then ! Particle in downstream
               ! find direction of v_2 from position x (radially outward):
               call radially_outward(phi_v, theta_v, x(1), x(2), x(3))

               ! v_2: velocity of downstream as seen in US
               v_2 = get_v_2(v_shock(t))
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               ! Lorentz transformed random step (advection)
               v_x = v_2*cos(phi_v)*sin(theta_v)
               v_y = v_2*sin(phi_v)*sin(theta_v)
               v_z = v_2*cos(theta_v)

               gamma_x = 1.d0/sqrt(1.d0 - v_x**2)                 ! v_2 dimless (v_2 = beta = v/c)
               gamma_y = 1.d0/sqrt(1.d0 - v_y**2)                 ! v_2 dimless (v_2 = beta = v/c)
               gamma_z = 1.d0/sqrt(1.d0 - v_z**2)                 ! v_2 dimless (v_2 = beta = v/c)

               x(1) = x(1) + v_x*dt + l_0*sin(theta)*cos(phi)/gamma_x
               x(2) = x(2) + v_y*dt + l_0*sin(theta)*sin(phi)/gamma_y
               x(3) = x(3) + v_z*dt + l_0*cos(theta)/gamma_z
            else ! Particle in upstream
               ! random step (isotropic in current medium rest frame)
               x(1) = x(1) + l_0*cos(phi)*sin(theta)
               x(2) = x(2) + l_0*sin(phi)*sin(theta)
               x(3) = x(3) + l_0*cos(theta)
            end if

            ! distances after step
            d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)              ! new particle distance
            t = t + dt
            E = E + dE
            r_sh2 = t_shock(t)                                  ! new shock dist

            ! Check for shock crossing
            if (d2 < r_sh2 .and. r_sh1 < d1) then ! Crossed shock (US->DS)
               call radially_outward(phi_v, theta_v, x(1), x(2), x(3)) ! angle of v_2 at pos x
               v_2 = get_v_2(v_shock(t)) ! US sees DS approach at velocity v_2
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               cos_theta = &
                  cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                  sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                  cos(theta)*cos(theta_v)

               gamma_v = 1.d0/sqrt(1.d0 - v_2**2)                 ! v_2 dimless (v_2 = beta = v/c)
               E_old = E                                          ! Old energy
               E = gamma_v*E*(1 - v_2*cos_theta)                  ! New energy
               rel_energy_gain = (E - E_old)/E_old                ! rel energy gain
               rel_energy_gain_sum = rel_energy_gain_sum + rel_energy_gain
               accel = 1
               num_crossings = num_crossings + 1
            else if (d2 > r_sh2 .and. r_sh1 > d1) then ! Cossed shock (DS -> US)
               call radially_outward(phi_v, theta_v, x(1), x(2), x(3)) ! direction of v_2 and shock
               v_2 = get_v_2(v_shock(t)) ! DS sees US approach at same velocity v_2
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               cos_theta = &
                  cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                  sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                  cos(theta)*cos(theta_v)

               gamma_v = 1.d0/sqrt(1.d0 - v_2**2)
               E_old = E
               E = gamma_v*E*(1 + v_2*cos_theta)
               rel_energy_gain = (E - E_old)/E_old
               rel_energy_gain_sum = rel_energy_gain_sum + rel_energy_gain
               accel = 1
               num_crossings = num_crossings + 1
            end if

            v_2 = get_v_2(v_shock(t))
            dmax = 3.d0*l_0_0/v_2
            f = f + df*dt                      ! \int dt f(t)
            delta = exp(-f)                    ! exp(-\int dt f(t))

            ! exit, if a) too late, b) too far down-stream, or c) scattering:
            if (t > t_max .or. d2 < r_sh2 - dmax .or. r > delta) then
               if (t > t_max .or. d2 < r_sh2 - dmax) then  ! we're tired or trapped behind
                  !              write(*,*) 'tired',n_in,n_out
                  r_sh0 = t_shock(t)
                  d0 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
                  if (t > t_max .and. d0 < r_sh0) then
                     print *, "Time exit - particle in downstream - num_cross: ", num_crossings
                  else if (t > t_max .and. d0 > r_sh0) then
                     print *, "Time exit - particle in upstream - num_cross: ", num_crossings
                  end if
                  call store(pid, E, w, num_crossings, rel_energy_gain_sum)
                  !call store_raw(E, set, n_injected, num_crossings)
                  print *, "Num shock crossings: ", num_crossings
                  print *, "Num steps taken: ", num_steps_taken
                  n_in = n_in - 1
                  n_out = n_out + 1
                  return
               end if
               if (r > delta) then                              ! decay or scattering
                  write (*, *) 'should never happen...'
                  stop
               end if
            end if

         end do
      end do
   end subroutine isotropic_random_walk

   subroutine scales_charged(m, Z, En, t, w, df, dt, dE)
      !use user_variables, only: t_max
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

   subroutine store(pid, En, w, num_crossings, rel_energy_gain_sum)
      use internal; use result, only: n_enbin, En_f, NE_esc, rel_energy_gain_total_sum
      implicit none
      integer pid, i
      double precision En, w, l
      integer, intent(in) :: num_crossings
      double precision, intent(in) :: rel_energy_gain_sum
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

      !rel_energy_gain_avg = rel_energy_gain_sum/num_crossings
      !rel_energy_gain_total_sum = rel_energy_gain_total_sum + rel_energy_gain_sum
      !  write(*,*) 'store: ',pid,i
   end subroutine store

   !subroutine store_raw(En, set_num, n_injected, num_crossings)
   !   ! Stores particles energies upon exit
   !   ! Raw, unbinned energies
   !   use result, only: exit_energies, num_crossings_total
   !   use user_variables, only: n_start
   !   implicit none
   !   double precision, intent(in) :: En
   !   integer, intent(in) :: set_num, n_injected, num_crossings
   !   integer :: idx
   !   !idx = n_injected + (set_num - 1)*n_start
   !   !exit_energies(idx) = En
   !   !num_crossings_total(idx) = num_crossings
   !end subroutine store_raw

end module acceleration
