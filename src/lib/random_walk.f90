! File: random_walk.f90

module random_walk
! Module containing the cone restricted random_walk procedure (future)
! 'scattering' angle restricted by theta_max (future)
! Shock acceleration can be turned on/off (future)
    implicit none
    private
    public pitch_angle_random_walk

contains
    subroutine pitch_angle_random_walk(set, n_injected) ! w/wo diffusion in trapping phase
        use user_variables, only: &
            debug, &
            theta_max, &
            shockless, &
            shockless_t_max, &
            num_sample_pos_target, &
            init_z, &
            no_stepsize_corr, &
            n_start
        use SNR_data, only: t_max; 
        use constants; use particle_data, only: m_p
        use event_internal; use result
        use internal
        use test_var, only: sec, accel

        implicit none
        integer, intent(in) :: set, n_injected
        integer :: k, n_step
        double precision :: r, m, f, df, dt, dE, delta, l_0, l_0_0
        double precision :: r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_rel
        double precision :: gamma_factor, cos_theta, dt_us
        double precision, dimension(3) :: l_vec, l_vec_us, v_rel_vec
        double precision :: ran0, R_L, t_shock, v_shock, get_v_rel
        integer, pointer :: pid, A, Z
        double precision, pointer :: E, x(:), t, w
        double precision, dimension(4, 4) :: boost_matrix
        double precision, dimension(4) :: l_four_vec
        double precision :: stepsize, analytical_stepsize, t0
        integer :: num_steps_taken, num_steps_total, sample_int, num_sample_pos, sample_count
        integer :: num_crossings

        pid => event(n_in)%pid
        A => event(n_in)%A
        Z => event(n_in)%Z
        E => event(n_in)%E
        x => event(n_in)%x
        t => event(n_in)%t
        w => event(n_in)%w

        d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
        if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6) then
            call error('wrong initial condition, shock', 0)
        end if
        r = ran0()
        m = A*m_p
        f = 0.d0

        if (shockless) then
            ! Initiate at (0, 0, 0)
            x(1) = 0.d0
            x(2) = 0.d0
            x(3) = 0.d0

            ! Find number of steps, sample interval etc. 
            t0 = t
            print *, "t0: ", t0
            if (no_stepsize_corr) then
                print *, "Isotropic stepsize (R_larmor)"
                l_0 =  R_L(E, t0)
            else
                l_0 = stepsize(E, t0, theta_max)
            end if
            dt = l_0
            if (l_0 <= 0.d0 .or. dt <= 0.d0) then
                print *, "Shockless: true"
                print *, "l_0: ", l_0
                call error("wrong stepsize", 0)
            end if
            num_steps_total = abs(t0 - shockless_t_max)/l_0 + 1
            sample_int = floor(real(num_steps_total/num_sample_pos_target)) + 1
            num_sample_pos = floor(real(num_steps_total/sample_int))
            sample_count = 0
            if (.not. allocated(sample_positions)) allocate(sample_positions(4, n_start, num_sample_pos))
        else
            num_crossings = 0
        end if
        num_steps_taken = 0

        do
            ! Trajectory and sample position log
            if (num_steps_taken + 1 <= size(initial_trajectories, 3)) then
                initial_trajectories(1, n_injected, num_steps_taken + 1) = t
                initial_trajectories(2, n_injected, num_steps_taken + 1) = x(1)
                initial_trajectories(3, n_injected, num_steps_taken + 1) = x(2)
                initial_trajectories(4, n_injected, num_steps_taken + 1) = x(3)
            end if
            if (shockless) then
                if (modulo(num_steps_taken + 1, sample_int) == 0) then
                    if (sample_count > num_sample_pos) then
                        print *, "sample_count: ", sample_count
                        print *, "num_sample_pos: ", num_sample_pos
                        call error("sample_count > num_sample_pos", 0)
                    end if
                    sample_count = sample_count + 1
                    sample_positions(1, n_injected, sample_count) = t
                    sample_positions(2, n_injected, sample_count) = x(1)
                    sample_positions(3, n_injected, sample_count) = x(2)
                    sample_positions(4, n_injected, sample_count) = x(3)
                end if
            end if

            if (.not. shockless) then
                ! Stepsize
                df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)            ! interaction rate (1/yr)
                call scales_charged(m, Z, E, t, w, df, dt, dE)
                !l_0 = R_L(E, t)/dble(Z)
                !l_0 = stepsize(E, t, theta_max)!R_L(E, t)/dble(Z)
                l_0 = analytical_stepsize(E, t, theta_max)!R_L(E, t)/dble(Z)
                l_0_0 = l_0
                if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)
            end if

            ! Step direction
            if (num_steps_taken == 0) then
                ! First step isotropic
                call isotropic_scatter(theta, phi)
            else
                ! Following steps are small angle
                call small_angle_scatter(theta, phi, theta_max)
            end if
            if (dt >= l_0) then                       ! one random step of size l_0
                dE = dE*l_0/dt
                dt = l_0
                n_step = 1
            else                                    ! n steps l0 in same direction
                ! Only relevant with interactions turned on
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
            num_steps_taken = num_steps_taken + 1
            !print *, x(1), x(2), x(3)

            ! Convert step from spherical to cartesian coordinates
            ! l_vec: cartesian step in the particle's local region's restframe
            ! l_0: stepsize in particle's local region's frame
            ! theta, phi: particle flight anlges in the local frame
            call spherical_to_cartesian(l_0, theta, phi, l_vec(1), l_vec(2), l_vec(3))

            do k = 1, n_step
                r_sh1 = t_shock(t)                              ! Radial position of shock before step
                d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! Radial position of particle before step

                if (d1 < r_sh1 .and. .not. shockless) then
                    ! Particle is in the downstream region
                    call cartesian_to_spherical(x(1), x(2), x(3), v_rel, theta_v, phi_v)
                    v_rel = get_v_rel(v_shock(t))
                    if (v_rel <= 0.d0) call error('v_rel<=0', 0)
                    call spherical_to_cartesian(v_rel, theta_v, phi_v, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))

                    ! Advection -- add downstream flow velocity to the step (Gallilean boost)
                    !l_vec(1) = l_vec(1) + v_rel_vec(1)*dt !v_rel*cos(phi_v)*sin(theta_v)*dt
                    !l_vec(2) = l_vec(2) + v_rel_vec(2)*dt !v_rel*sin(phi_v)*sin(theta_v)*dt
                    !l_vec(3) = l_vec(3) + v_rel_vec(3)*dt !v_rel*cos(theta_v)*dt

                    ! Advection -- account for downstream velocity -- Lorentz transform
                    call lorentz_boost_matrix(-v_rel_vec, boost_matrix)
                    l_four_vec = matmul(boost_matrix, [dt, l_vec(1), l_vec(2), l_vec(3)])
                    dt_us = l_four_vec(1)
                    l_vec_us = l_four_vec(2:)

                    dt = dt_us
                    l_vec(1) = l_vec_us(1)
                    l_vec(2) = l_vec_us(2)
                    l_vec(3) = l_vec_us(3)
                end if
                x(1) = x(1) + l_vec(1)
                x(2) = x(2) + l_vec(2)
                x(3) = x(3) + l_vec(3)

                t = t + dt
                d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)              ! new distance
                E = E + dE
                r_sh2 = t_shock(t)

                if (.not. shockless) then
                if (d2 < r_sh2 .and. d1 > r_sh1) then
                    ! we have crossed to the left: US -> DS
                    call cartesian_to_spherical(x(1), x(2), x(3), v_rel, theta_v, phi_v)
                    v_rel = get_v_rel(v_shock(t))
                    if (v_rel <= 0.d0) call error("v_rel <= 0", 0)
                    gamma_factor = 1.d0/(1.d0 - v_rel**2)
                    cos_theta = &
                        cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                        sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                        cos(theta)*cos(theta_v)
                    E = gamma_factor*E*(1.d0 - v_rel*cos_theta)
                    accel = 1
                else if (d2 > r_sh2 .and. d1 < r_sh1) then
                    ! we have crossed to the right: DS -> US
                    call cartesian_to_spherical(x(1), x(2), x(3), v_rel, theta_v, phi_v)
                    v_rel = get_v_rel(v_shock(t))
                    if (v_rel <= 0.d0) call error("v_rel <= 0", 0)
                    gamma_factor = 1.d0/(1.d0 - v_rel**2)
                    cos_theta = &
                        cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                        sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                        cos(theta)*cos(theta_v)
                    E = gamma_factor*E*(1.d0 + v_rel*cos_theta)
                    accel = 1
                end if
                end if

                v_rel = 0.75d0*v_shock(t)
                dmax = 3.d0*l_0_0/v_rel
                f = f + df*dt                      ! \int dt f(t)
                delta = exp(-f)                    ! exp(-\int dt f(t))

                ! shockless exit
                if (shockless) then
                    if (t > shockless_t_max) then
                        ! exit random walking particle when max time exceeded
                        print *, "shockless_t_max: ", shockless_t_max
                        print *, "t: ", t
                        call store_shockless(n_injected, x(1), x(2), x(3))
                        print *, "t-max: ", shockless_t_max
                        print *, "num steps tot: ", num_steps_total
                        print *, "num steps taken: ", num_steps_taken
                        n_in = n_in - 1
                        n_out = n_out + 1
                        return
                    end if
                ! acceleration exit
                else if (t > t_max .or. d2 < r_sh2 - dmax .or. r > delta) then
                    ! exit accel particle, if a) too late, b) too far down-stream, or c) scattering:
                    if (t > t_max .or. d2 < r_sh2 - dmax) then  ! we're tired or trapped behind
                        !              write(*,*) 'tired',n_in,n_out
                        call store(pid, E, w)
                        !print *, "Exit energy: ", E
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
    end subroutine pitch_angle_random_walk

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
        !tausyn = tau_syn(m,En,t)/dble(Z**2)                           ! synchrotron
        tausyn = huge(0.d0)
        dt = 9.d-3*min(tau_eff, tausyn)

        dE_loss_syn = dt/tausyn*En                                             ! eV
        dE = -dE_loss_syn

        if (abs(dE)/En > 1.d-2) call error('dE too large in scale', 1)
        ! no sync. photo-production
    end subroutine scales_charged

    subroutine store(pid, En, w)
        use internal; use result

        implicit none
        integer pid, i
        double precision En, w, l

        l = log10(En)                                                ! energy bin
        i = int((l - d_f)/dn)
        if (i <= 0) then
            call error('stor_esc, E<=Emin', 11)
            i = 1
        end if
        if (i > n_enbin) then
            i = n_enbin
        end if
        En_f(pid, i) = En_f(pid, i) + w*En
    end subroutine store

    subroutine store_shockless(n_injected, x1, x2, x3)
        use result, only: final_positions
        implicit none
        double precision, intent(in) :: x1, x2, x3
        integer, intent(in) :: n_injected
        
        final_positions(1, n_injected) = x1
        final_positions(2, n_injected) = x2
        final_positions(3, n_injected) = x3
    end subroutine store_shockless

end module random_walk
