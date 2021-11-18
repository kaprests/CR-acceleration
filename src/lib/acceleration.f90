module non_rel_acceleration
   implicit none
   private

   public isotropic_random_walk
contains


   subroutine isotropic_random_walk(set, n_injected)
      ! Isotropic random walk
      use user_variables, only: debug, t_max
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
      double precision :: gamma_v, cos_theta
      double precision :: d0, r_sh0

      ! ############
      integer :: num_crossings
      double precision :: rel_energy_gain, E_old, rel_energy_gain_sum

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

      rel_energy_gain_sum = 0
      num_crossings = 0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Simulate DSA until exit !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do
         !!! Step size, step number etc. !!!

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

         ! find step size and new position:
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

         !!! perform random step(s) !!!
         do k = 1, n_step
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!! Diffusion/advection !!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! distances before step
            r_sh1 = t_shock(t)                              ! old shock position
            d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! old distance/particle radial position

            ! Perform random step
            if (d1 < r_sh1) then ! Particle in downstream
               ! find direction of v_2 from position x (radially outward):
               call radially_outward(phi_v,theta_v, x(1), x(2), x(3))

               ! v_2: velocity of downstream as seen in US
               v_2 = get_v_2(v_shock(t))
               if (v_2 <= 0.d0) call error('v_2<=0', 0)

               ! Lorentz transformed random step (advection)
               gamma_v = 1.d0/sqrt(1.d0 - v_2**2)                 ! v_2 dimless (v_2 = beta = v/c)
               x(1) = x(1) + v_2*cos(phi_v)*sin(theta_v)*dt + l_0*sin(theta)*cos(phi)/gamma_v
               x(2) = x(2) + v_2*sin(phi_v)*sin(theta_v)*dt + l_0*sin(theta)*sin(phi)/gamma_v
               x(3) = x(3) + v_2*cos(theta_v)*dt + l_0*cos(theta)/gamma_v
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

            !!!!!!!!!!!!!!!!!!!!!!
            !!! Shock crossing !!!
            !!!!!!!!!!!!!!!!!!!!!!
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
                  d0 = sqrt(x(1)**2+x(2)**2+x(3)**2)
                  if (t > t_max .and. d0 < r_sh0) then
                     print *, "Time exit - particle in downstream - num_cross: ", num_crossings
                  else if (t > t_max .and. d0 > r_sh0) then
                     print *, "Time exit - particle in upstream - num_cross: ", num_crossings
                  end if
                  call store(pid, E, w, num_crossings, rel_energy_gain_sum)
                  call store_raw(E, set, n_injected)
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


   subroutine store(pid, En, w, num_crossings, rel_energy_gain_sum)
      use internal; use result, only : n_enbin, En_f, NE_esc, rel_energy_gain_total_sum
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
      rel_energy_gain_avg = rel_energy_gain_sum / num_crossings
      rel_energy_gain_total_sum = rel_energy_gain_total_sum + rel_energy_gain_sum
   !  write(*,*) 'store: ',pid,i
   end subroutine store


   subroutine store_raw(En, set_num, particle_num)
      ! Stores particles energies upon exit
      ! Raw, unbinned energies
      use result
      implicit none
      double precision, intent(in) :: En
      integer, intent(in) :: set_num, particle_num
      exit_energies(set_num, particle_num) = En
   end subroutine store_raw


   subroutine anisotropic_upstream(phi, theta, x1, x2, x3)
      ! Angles of anisotropized upstream particles
      ! Limited, random deflection centered around the shock normal
      use constants, only : pi, two_pi
      use user_variables, only: gamma_sh
      implicit none
      double precision, intent(in) :: x1, x2, x3
      double precision, intent(inout) :: phi, theta
      double precision :: max_deflection, theta_def, phi_def
      double precision :: ran0
      max_deflection = 2/gamma_sh
      call radially_outward(phi, theta, x1, x2, x3)
      theta_def = -max_deflection + 2*max_deflection*ran0()
      phi_def = -max_deflection + 2*max_deflection*ran0()
      theta = theta + theta_def
      phi = phi + phi_def
      if (theta < 0) then
         theta = -theta
         phi = phi + pi
      else if (theta > pi) then
         theta = two_pi - theta
         phi = phi + pi
      endif
      phi = modulo(phi, two_pi)
   end subroutine anisotropic_upstream


end module non_rel_acceleration
