module random_walk
    implicit none
    private

    public pitch_angle_random_walk !, isotropic_random_walk
contains

    subroutine pitch_angle_random_walk(set, n_injected)
        ! Pitch angle scattering shock acceleration
        use user_variables, only: &
            debug, &
            t_max, &
            num_steps_log, &
            inj_model, &
            shockless, &
            num_steps_log, &
            num_phase_log, &
            z_axis, &
            no_small_angle_corr
        use constants
        use particle_data, only: m_p
        use event_internal
        use result
        use internal
        use test_var, only: sec, accel

        implicit none
        integer, intent(in) :: set, n_injected          ! Dummy variables
        integer :: k, n_step                            ! Loop variables: incr, end
        double precision :: &                           ! Functions
            ran0, &
            R_L, &
            t_shock, &
            v_shock, &
            get_v_2, &
            stepsize, &
            analytical_stepsize, &
            spacetime_interval
        integer, pointer :: pid, A, Z                   ! Integer pointers
        double precision, pointer :: E, x(:), t, w      ! DP pointers
        double precision :: &
            r_sh1, r_sh2, &                              ! Radial pos of shockfront bef/aft step
            phi, theta, &                                ! Random step angles (direction)
            phi_v, theta_v, &                            ! Shock normal angles
            d1, d2, dmax, &                              ! Radial pos of particle bef/aft step
            v_2, v_p, &                                  ! Relative velocity US/DS shock and particle
            dist_particle_shock, &                       ! Distance between particle and shock
            p_particle, &
            E_srf                                        ! Energy in shock rest frame
        double precision :: &
            r, m, f, df, dt, dE, delta, &
            l_0, l_0_0, gamma_v
        double precision :: v_x, v_y, v_z               ! Cartesian comp of shock velocity
        double precision :: v_px, v_py, v_pz            ! Cartesian comp of particle velocity
        double precision :: px, py, pz                  ! Cartesian comp of particle 3-momentum
        double precision :: &
            cos_theta                                    ! Flight angle cosine
        double precision :: d0, r_sh0                   ! Initial particle and shock radial positions
        integer :: num_crossings, num_steps_taken
        double precision :: &
            rel_energy_gain, E_old, rel_energy_gain_sum
        double precision :: theta_max                   ! Max scattering angles
        double precision :: &
            g(3), p(3), R_euler(3, 3), theta_e, phi_e    ! Pitch angle quantities
        integer :: i, j, idx, l                          ! Iteration variables
        double precision :: t0
        integer :: sample_int, sample_count, num_samples, num_steps_total
        double precision, dimension(3) :: l_vec        ! step vector -- cartesian coordinates
        double precision :: dt_us                       ! upstream time step for a downstream particle
        double precision, dimension(3) :: l_vec_us      ! upstream spatial step for a downstream particle

        pid => event(n_in)%pid
        A => event(n_in)%A
        Z => event(n_in)%Z
        E => event(n_in)%E
        x => event(n_in)%x
        t => event(n_in)%t
        w => event(n_in)%w

        ! Initial radial position of particle
        if (.not. shockless) then
            d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
            if (sec == 0 .and. abs(d1/t_shock(t) - 1.01d0) .gt. 1.d-6 .and. inj_model == 0) then
                call error('wrong initial condition, shock', 0)
            else if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6 .and. inj_model == 1) then
                call error('wrong initial condition, shock', 0)
            else if (sec == 0 .and. abs(d1/t_shock(t) - 1.d0) .gt. 1.d-6 .and. inj_model == 2) then
                call error('wrong initial condition, shock', 0)
            end if
        end if

        if (t < t_0_0) then
            print *, "t: ", t
            print *, "t_0_0: ", t_0_0
            call error("t < t_0_0 (initial) -- wrong t_0_0", 0)
        end if

        r = ran0()
        m = A*m_p
        f = 0.d0

        ! Set theta_max (max scattering angle) for whole simulation run
        call set_theta_max(theta_max)

        if (shockless) then
            ! Initiate at (0, 0, 0)
            x(1) = 0.d0
            x(2) = 0.d0
            x(3) = 0.d0
            ! Find number of steps, sample intervals etc.
            t0 = t
            if (no_small_angle_corr) then
                print *, "ISO"
                l_0 = R_L(E, t0)
            else
                l_0 = stepsize(E, t0, theta_max)
            end if
            dt = l_0
            if (l_0 <= 0.d0 .or. dt <= 0.d0) then
                write(*, *) "SHOCKLESS"
                write(*, *) "l_0: ", l_0
                call error('wrong scales', 0)
            end if
            num_steps_total = abs(t0 - t_max)/l_0 + 1
            sample_int = floor(real(num_steps_total/num_steps_log)) + 1 ! sample interval
            num_samples = floor(real(num_steps_total/sample_int))
            sample_count = 0
            if (.not. allocated(sample_positions)) allocate (sample_positions(4, n_start, num_samples))
        else
            rel_energy_gain_sum = 0
            num_crossings = 0
        end if

        ! Both shock on and off
        num_steps_taken = 0

        do
            t0 = t ! Time before step
            !!!!!!!!!!!!!!!!
            ! Log position ! (initial trajectories)
            !!!!!!!!!!!!!!!!
            if (num_steps_taken + 1 <= size(trajectories, 3)) then
                trajectories(1, n_injected, num_steps_taken + 1) = x(1)
                trajectories(2, n_injected, num_steps_taken + 1) = x(2)
                trajectories(3, n_injected, num_steps_taken + 1) = x(3)
                trajectories(4, n_injected, num_steps_taken + 1) = t
            end if

            !!!!!!!!!!!!!!!!!!!!
            ! Sample positions !
            !!!!!!!!!!!!!!!!!!!!
            if (shockless) then
                if (modulo(num_steps_taken + 1, sample_int) == 0) then
                    if (sample_count > num_samples) then
                        write (*, *) "sample_count: ", sample_count
                        write (*, *) "num_samples: ", num_samples
                        write (*, *) "num_int: ", sample_int
                        write (*, *) "num_steps_taken: ", num_steps_taken
                        write (*, *) "num_steps_total: ", num_steps_total
                        call error('sample count too big!', 1)
                    end if
                    sample_count = sample_count + 1
                    sample_positions(1, n_injected, sample_count) = x(1)
                    sample_positions(2, n_injected, sample_count) = x(2)
                    sample_positions(3, n_injected, sample_count) = x(3)
                    sample_positions(4, n_injected, sample_count) = t
                end if
            end if

            !!!!!!!!!!!!
            ! Stepsize !
            !!!!!!!!!!!!
            df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)   ! interaction rate (1/yr)
            call scales_charged(m, Z, E, t, w, df, dt, dE)
            if (no_small_angle_corr) then
                ! Isotropic rw stepsize
                l_0 = R_L(E, t)
            else
                ! Small angle stepsize 
                l_0 = stepsize(E, t, theta_max)
            end if
            l_0_0 = l_0
            if (l_0 <= 0.d0 .or. dt <= 0.d0) then
                write(*, *) "SHOCKLESS: ", shockless
                write(*, *) "l_0: ", l_0
                call error('wrong scales', 0)
            end if
            ! Number of steps
            if (dt >= l_0) then                         ! one random step of size l_0
                dE = dE*l_0/dt
                dt = l_0
                n_step = 1
            else                                        ! n steps l0 in same direction
                l_0 = dt
                if (l_0_0/dt < 1.d3) then
                    n_step = int(l_0_0/dt + 0.5d0)
                else                                    ! fast decays lead to overflow
                    n_step = 1000                       ! this should be enough
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

            !!!!!!!!!!!!!!!!!!
            ! Step direction !
            !!!!!!!!!!!!!!!!!!
            r_sh0 = t_shock(t)                              ! shock position 
            d0 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! particle radial position
            if (num_steps_taken == 0) then
                if (z_axis .and. shockless) then
                    ! First step along z-axis
                    theta = 0.0
                    phi = 0.0 ! Can be anything
                else
                    ! First step isotropic, always the case when there is a shock
                    call isotropic(theta, phi)
                end if
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
                        ! Todo: define R_euler as transposed (fortran is column major)
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

            ! Step vector l_vec
            call spherical_to_cartesian(l_0, theta, phi, l_vec(1), l_vec(2), l_vec(3))

            ! Increment number of steps taken
            num_steps_taken = num_steps_taken + 1

            !!!!!!!!!!!!!!!!!!!
            ! Perform step(s) !
            !!!!!!!!!!!!!!!!!!!
            do k = 1, n_step
                r_sh1 = t_shock(t)                          ! old shock position
                d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)      ! old distance/particle radial position

                ! Perform random step
                if (d1 < r_sh1 .and. .not. shockless) then 
                    ! Particle in downstream (advection)
                    ! => generated l_0, theta, phi and dt are downstream quantities
                    ! Need to tranform (boost) to labframe (upstream)

                    ! Boost the step vector to the lab (US) frame
                    call downstream_to_upstream_boost(dt, l_vec, dt_us, l_vec_us, v_shock(t), t, [x(1), x(2), x(3)])

                    ! Set the upstream values
                    dt = dt_us
                    l_vec = l_vec_us
                    l_0 = norm2(l_vec)
                    !call radially_outward(theta, phi, l_vec(1), l_vec(2), l_vec(3))
                end if
                t = t + dt
                x(1) = x(1) + l_vec(1)
                x(2) = x(2) + l_vec(2)
                x(3) = x(3) + l_vec(3)

                ! distances after step
                d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)              ! new particle distance
                E = E + dE
                r_sh2 = t_shock(t)                                  ! new shock dist

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Check for shock crossing !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (d2 < r_sh2 .and. r_sh1 < d1 .and. .not. shockless) then ! Crossed (US->DS)
                    call radially_outward(phi_v, theta_v, x(1), x(2), x(3)) ! angle of v_2 at pos x
                    v_2 = get_v_2(v_shock(t)) ! US sees DS approach at velocity v_2
                    if (v_2 <= 0.d0) call error('v_2<=0', 0)

                    cos_theta = &
                        cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &
                        sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &
                        cos(theta)*cos(theta_v)

                    gamma_v = 1.d0/sqrt(1.d0 - v_2**2)            ! v_2 dimless (v_2 = beta = v/c)
                    E_old = E                                     ! Old energy
                    E = gamma_v*E*(1 - v_2*cos_theta)             ! New energy
                    rel_energy_gain = (E - E_old)/E_old           ! rel energy gain
                    rel_energy_gain_sum = rel_energy_gain_sum + rel_energy_gain
                    accel = 1
                    num_crossings = num_crossings + 1

                    ! log angles at crossing
                    !if (num_crossings <= size(crossing_flight_angles, 2)) then
                    !    crossing_flight_angles(n_injected, num_crossings) = acos(cos_theta)
                    !end if

                    ! Lorentz tranform three momentum
                else if (d2 > r_sh2 .and. r_sh1 > d1 .and. .not. shockless) then ! Cossed (DS -> US)
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

                    !if (num_crossings <= size(crossing_flight_angles, 2)) then
                    !    crossing_flight_angles(n_injected, num_crossings) = acos(cos_theta)
                    !end if

                    ! Lorentz tranform three momentum
                end if

                v_2 = get_v_2(v_shock(t))
                dmax = 3.d0*l_0_0/v_2
                f = f + df*dt                      ! \int dt f(t)
                delta = exp(-f)                    ! exp(-\int dt f(t))

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Exit -- Random walking particle !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (shockless .and. t > t_max) then
                    idx = n_injected + (set - 1)*n_start
                    if (idx == n_start*n_sets) then
                        !print *, "Particle exiting"
                        !print *, "t max: ", t_max
                        !print *, "steps taken: ", num_steps_taken
                        !print *, "num_steps_total: ", num_steps_total
                        !print *, "num_steps_log: ", num_steps_log
                        !print *, "sample_count: ", sample_count
                        !print *, "num_samples: ", num_samples
                    end if
                    call store_shockless(n_injected, x(1), x(2), x(3))
                    n_in = n_in - 1
                    n_out = n_out + 1
                    return
                end if

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Exit -- accelerating particle !
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if (.not. shockless) then
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
                            !print *, "!!!!!!!!!!!!!!!!!"
                            !print *, "w: ", w
                            !print *, "!!!!!!!!!!!!!!!!!"
                            call store(pid, E, w, num_crossings, rel_energy_gain_sum)
                            !!call store_raw(E, set, n_injected, num_crossings)
                            !print *, "#############################"
                            print *, "Num shock crossings before exit: ", num_crossings
                            !print *, "Num steps taken: ", num_steps_taken
                            !print *, "Final theta max: ", theta_max
                            !print *, "t-exit: ", t
                            !print *, "t-max: ", t_max
                            !print *, "t_max/100*5: ", t_max/500
                            n_in = n_in - 1
                            n_out = n_out + 1
                            return
                        end if
                        if (r > delta) then                              ! decay or scattering
                            write (*, *) 'should never happen...'
                            stop
                        end if
                    end if
                end if
            end do
        end do
    end subroutine pitch_angle_random_walk

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
        ! tausyn = tau_syn(m,En,t)/dble(Z**2)                           ! synchrotron
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
        !write(*,*) 'store: ',pid,i
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
