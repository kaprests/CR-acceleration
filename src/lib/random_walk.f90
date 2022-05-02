! File: random_walk.f90

module random_walk
! Module containing the cone restricted random_walk procedure 
! 'scattering' angle restricted by theta_max 
! Shock acceleration can be turned on/off 
    implicit none
    private
    public pitch_angle_random_walk

contains
    subroutine pitch_angle_random_walk(set, n_injected, n_proc) ! w/wo diffusion in trapping phase
        use user_variables, only: &
            debug, &
            theta_max, &
            shockless, &
            shockless_t_max, &
            num_sample_pos_target, &
            init_z, &
            no_stepsize_corr, &
            n_start, &
            n_sets, &
            inj_model, &
            gamma_shock_const
        use SNR_data, only: t_max; 
        use constants; use particle_data, only: m_p
        use event_internal; use result
        use internal
        use test_var, only: sec, accel
        use result

        implicit none
        integer, intent(in) :: set, n_injected, n_proc
        integer :: k, n_step
        double precision :: r, m, f, df, dE, delta, l_0, l_0_0
        double precision :: r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_rel
        double precision :: gamma_factor, gamma_particle, cos_theta
        double precision :: dt_loc, dt_us, dt_ds
        double precision, dimension(3) :: l_vec_loc, l_vec_us, l_vec_ds, v_rel_vec
        double precision :: ran0, R_L, t_shock, v_shock, get_v_rel
        integer, pointer :: pid, A, Z
        double precision, pointer :: E, x(:), t, w, p(:)
        double precision :: stepsize, analytical_stepsize, t0, v_particle
        double precision :: upstream_loss_cone_opening_angle
        integer :: num_steps_taken, num_steps_total, sample_int, num_sample_pos, sample_count
        integer :: num_crossings
        double precision :: E_0, v_p, cross_angle
        double precision, dimension(3) :: p_0

        pid => event(n_in)%pid
        A => event(n_in)%A
        Z => event(n_in)%Z
        E => event(n_in)%E
        x => event(n_in)%x
        p => event(n_in)%p
        t => event(n_in)%t
        w => event(n_in)%w

        d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
        if (sec == 0 .and. abs(d1/(t_shock(t)+stepsize(E, t, theta_max)) - 1.d0) .gt. 1.d-6) then
            call error('wrong initial condition, shock', 0)
        end if
        r = ran0()
        m = A*m_p
        f = 0.d0

        if (shockless) then
            ! Initiate at origin
            x(1) = 0.d0
            x(2) = 0.d0
            x(3) = 0.d0

            ! Find number of steps, sample interval etc. 
            t0 = t
            if (no_stepsize_corr) then
                l_0 =  R_L(E, t0)
            else
                l_0 = stepsize(E, t0, theta_max)
            end if
            dt_loc = l_0/v_particle(E, m_p)
            if (l_0 <= 0.d0 .or. dt_loc <= 0.d0) then
                print *, "Shockless: true"
                print *, "l_0: ", l_0
                call error("wrong stepsize", 0)
            end if
            num_steps_total = abs(t0 - shockless_t_max)/l_0 + 1
            sample_int = floor(real(num_steps_total/num_sample_pos_target)) + 1
            num_sample_pos = floor(real(num_steps_total/sample_int))
            sample_count = 0
            if (.not. allocated(sample_positions)) &
                allocate(sample_positions(4, n_start, num_sample_pos))
        else
            num_crossings = 0
        end if
        num_steps_taken = 0

        do
            ! Log position for initial trajectory and sampled positions
            if (num_steps_taken + 1 <= size(initial_trajectories, 3)) then
                ! Initial trajectory
                initial_trajectories(1, n_injected, num_steps_taken + 1) = t
                initial_trajectories(2, n_injected, num_steps_taken + 1) = x(1)
                initial_trajectories(3, n_injected, num_steps_taken + 1) = x(2)
                initial_trajectories(4, n_injected, num_steps_taken + 1) = x(3)
            end if
            if (shockless) then
                if (modulo(num_steps_taken + 1, sample_int) == 0) then
                    ! Sampled position
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
                ! Update stepsize for accelerating particle
                df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)            ! interaction rate (1/yr)
                call scales_charged(m, Z, E, t, w, df, dt_loc, dE)
                l_0 = stepsize(E, t, theta_max)!R_L(E, t)/dble(Z)
                l_0_0 = l_0
                dt_loc = l_0/v_particle(E, m_p)
                if (l_0 <= 0.d0 .or. dt_loc <= 0.d0) call error('wrong scales', 0)
            end if

            ! Step direction
            if (num_steps_taken == 0) then
                ! First step isotropic
                call isotropic_scatter(theta, phi)
                if (init_z) then
                    ! First step along z-axis
                    theta = 0.d0
                end if
            else
                ! Following steps are small angle
                call small_angle_scatter(theta, phi, theta_max)
            end if

            if (dt_loc >= l_0) then                       ! one random step of size l_0
                dE = dE*l_0/dt_loc
                n_step = 1
            else                                    ! n steps l0 in same direction
                ! Only relevant with interactions turned on
                l_0 = dt_loc * v_particle(E, m_p)
                if (l_0_0/dt_loc < 1.d3) then
                    n_step = int(l_0_0/dt_loc + 0.5d0)
                else                                 ! fast decays lead to overflow
                    n_step = 1000                    ! this should be enough
                end if
                if (debug > 0) write (*, *) 'E, step number', E, n_step
            end if
            if (n_step < 1) then
                write (*, *) l_0_0, l_0
                write (*, *) dt_loc, df
                write (*, *) A, Z
                write (*, *) E
                write (*, *) l_0_0/dt_loc, n_step
                call error('wrong step number', 0)
            end if
            num_steps_taken = num_steps_taken + 1
            if (n_sets * n_start * n_proc == 1) then
                print *, x(1), x(2), x(3)
            end if

            ! Step and momentum in the particle's local region rest frame

            ! l_vec: spatial step vector in local frame
            call spherical_to_cartesian(&
                l_0, theta, phi, l_vec_loc(1), l_vec_loc(2), l_vec_loc(3))

            ! gamma factor of the particle
            gamma_particle = 1.d0/sqrt(1.d0 - v_particle(E, m_p)**2)

            ! p: three-momentum vector for particle in the local frame
            call spherical_to_cartesian(&
                gamma_particle*m_p* &
                sqrt(l_vec_loc(1)**2 + l_vec_loc(2)**2 + l_vec_loc(3)**2)/dt_loc, &
                theta, phi, p(1), p(2), p(3))

            ! if particle in downstream -> transform step to upstream
            ! Then apply step, propagate time by dt (upstream), and find new positions
            do k = 1, n_step
                r_sh1 = t_shock(t)  ! Radial position of shock before step 
                d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)  ! Radial position of particle before step

                ! transform step if particle is in the downstream region
                if (.not. shockless) then
                    if (d1 < r_sh1) then ! Particle is in the downstream region
                        dt_ds = dt_loc          ! downstream timestep
                        l_vec_ds = l_vec_loc    ! downstream step vector
                        call cartesian_to_spherical(& ! Radial out direction
                            x(1), x(2), x(3), v_rel, theta_v, phi_v)
                        v_rel = get_v_rel(v_shock(t)) ! Relative velocity of DS and US
                        if (v_rel <= 0.d0) call error('v_rel<=0', 0)
                        ! relative 3-velocity vector:
                        call spherical_to_cartesian(& ! shock 3-velocity in downstream frame
                            v_rel, theta_v, phi_v, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
                        v_rel_vec = -v_rel_vec ! Upstream flows radially inward (as seen in DS)
                        ! dt_us and l_vec_us: 
                        call lorentz_boost(dt_ds, l_vec_ds, dt_us, l_vec_us, v_rel_vec)
                    else
                        dt_us = dt_loc          ! dt_us
                        l_vec_us = l_vec_loc    ! l_vec_us
                        ! Don't need dt_ds and l_vec_ds in this case
                    end if
                end if
                t = t + dt_us
                x = x + l_vec_us

                ! Update particle distance and shock positions
                r_sh2 = t_shock(t)
                d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)              ! new distance
                E = E + dE

                if (.not. shockless) then
                if (d2 < r_sh2 .and. d1 >= r_sh1) then ! crossed to the 'left': US -> DS
                    ! Radial (out) direction:
                    call cartesian_to_spherical(& 
                        x(1), x(2), x(3), v_rel, theta_v, phi_v)

                    ! Relative velocity of DS and US:
                    v_rel = get_v_rel(v_shock(t))
                    if (v_rel <= 0.d0) call error('v_rel<=0', 0)

                    ! relative 3-velocity vector:
                    call spherical_to_cartesian(& ! shock 3-velocity in downstream frame
                        v_rel, theta_v, phi_v, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))

                    E_0 = E ! local frame before crossing i.e. upstream frame
                    p_0 = p ! local frame before crossing i.e. upstream frame

                    ! E, p: in local frame after crossing i.e. downstream:
                    call lorentz_boost(E_0, p_0, E, p, v_rel_vec)

                    ! Store cross angle -- angle measured in upstream frame:
                    cross_angle = & ! planar shock approx, more accurate for smaller steps
                        acos(dot_product(l_vec_us, [x(1), x(2), x(3)]) / (l_0 * d1))
                    if (num_crossings == 0) then
                        call store_cross_angle_iso(cross_angle, cross_angle_distribution_first)
                    else
                        !call store_angle(cross_angle, cross_angle_distribution_updown, pi)
                        call store_cross_angle_iso(cross_angle, cross_angle_distribution_updown)
                        if (inj_model == 0 .and. gamma_shock_const >= 100.d0) then
                            call store_cross_angle_aniso(&
                                cross_angle, &
                                cross_angle_distribution_aniso_updown &
                            )
                        end if
                    end if

                    ! flight direction in new frame:
                    call cartesian_to_spherical(p(1), p(2), p(3), v_p, theta, phi)

                    accel = 1
                    num_crossings = num_crossings + 1
                else if (d2 >= r_sh2 .and. d1 < r_sh1) then ! crossed to the right: DS -> US
                    ! Radial (out) direction:
                    call cartesian_to_spherical(& 
                        x(1), x(2), x(3), v_rel, theta_v, phi_v)

                    ! Relative velocity of DS and US:
                    v_rel = get_v_rel(v_shock(t))
                    if (v_rel <= 0.d0) call error('v_rel<=0', 0)

                    ! relative 3-velocity vector:
                    call spherical_to_cartesian(& ! shock 3-velocity in downstream frame
                        v_rel, theta_v, phi_v, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
                    v_rel_vec = -v_rel_vec ! US flows radially inward

                    E_0 = E ! E_0 energy in local frame before cross i.e. downstream frame
                    p_0 = p ! p_0 energy in local frame before cross i.e. downstream frame
                    ! E, p: energy and momentum in local frame after cross i.e. upstream
                    call lorentz_boost(E_0, p_0, E, p, v_rel_vec)

                    ! Store cross angle -- angle measured in DS frame
                    cross_angle = & ! approx, more accurate for smaller steps
                        acos(dot_product(l_vec_ds, [x(1), x(2), x(3)])/(l_0 * d1))
                    call store_cross_angle_iso(cross_angle, cross_angle_distribution_downup)

                    ! flight direction in new frame: 
                    call cartesian_to_spherical(p(1), p(2), p(3), v_p, theta, phi)

                    accel = 1
                    num_crossings = num_crossings + 1
                end if
                end if

                v_rel = get_v_rel(v_shock(t))
                dmax = 3.d0*R_L(E, t)/v_rel

                ! Interaction stuff: 
                ! NB: perhaps dt_us not dt_loc in assignment of f?
                f = f + df*dt_loc                  ! \int dt f(t) 
                delta = exp(-f)                    ! exp(-\int dt f(t))

                ! shockless exit:
                if (shockless) then
                    if (t > shockless_t_max) then
                        ! exit random walking particle when max time exceeded
                        call store_shockless(n_injected, x(1), x(2), x(3))
                        n_in = n_in - 1
                        n_out = n_out + 1
                        return
                    end if
                ! acceleration exit:
                else if (t > t_max .or. d2 < r_sh2 - dmax .or. r > delta) then
                    ! exit accel particle, if a) too late, b) too far down-stream, or c) scattering:
                    if (t > t_max .or. d2 < r_sh2 - dmax) then  ! we're tired or trapped behind
                        !              write(*,*) 'tired',n_in,n_out
                        if (d2 < r_sh2 - dmax) then
                            ! particle exited in DS - transform energy to US
                            ! Radial (out) direction:
                            call cartesian_to_spherical(& 
                                x(1), x(2), x(3), v_rel, theta_v, phi_v)

                            ! Relative velocity of DS and US:
                            v_rel = get_v_rel(v_shock(t))
                            if (v_rel <= 0.d0) call error('v_rel<=0', 0)

                            ! relative 3-velocity vector:
                            call spherical_to_cartesian(& ! shock 3-velocity in downstream frame
                                v_rel, theta_v, phi_v, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
                            v_rel_vec = -v_rel_vec ! US flows radially inward

                            E_0 = E ! E_0 energy in local frame before cross i.e. downstream frame
                            p_0 = p ! p_0 energy in local frame before cross i.e. downstream frame
                            ! E, p: energy and momentum in local frame after cross i.e. upstream
                            call lorentz_boost(E_0, p_0, E, p, v_rel_vec)
                            print *, "DS exit, num crossings: ", num_crossings
                        else
                            print *, "Upstream time exit, num crossings: ", num_crossings
                        end if
                        if (n_start * n_sets * n_proc == 1) then
                            print *, "num crossings: ", num_crossings
                        end if
                        call store(pid, E, w)
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
        exit_energy_enumerate_dist(i) = exit_energy_enumerate_dist(i) + 1
    end subroutine store

    subroutine store_cross_angle_iso(angle, distribution_array)
        use constants, only: pi
        use result, only: n_angle_bins
        implicit none
        double precision, intent(in) :: angle
        double precision, intent(inout) :: distribution_array(n_angle_bins)

        call store_angle(angle, distribution_array, n_angle_bins, pi)
    end subroutine store_cross_angle_iso

    subroutine store_cross_angle_aniso(angle, distribution_array)
        use constants, only: pi
        use result, only: n_angle_bins_aniso
        implicit none
        double precision, intent(in) :: angle
        double precision, intent(inout) :: distribution_array(n_angle_bins_aniso)

        call store_angle(angle, distribution_array, n_angle_bins_aniso, pi/2.d0)
    end subroutine store_cross_angle_aniso

    subroutine store_angle(angle, distribution_array, n_angle_bins, max_angle)
        use internal

        implicit none
        double precision, intent(in) :: angle, max_angle
        integer, intent(in) :: n_angle_bins
        double precision, intent(inout) :: distribution_array(n_angle_bins)
        double precision :: bin_size
        integer :: bin_num

        bin_size = max_angle/n_angle_bins
        bin_num = max(ceiling(angle/bin_size), 1)
        distribution_array(bin_num) = distribution_array(bin_num) + 1
    end subroutine store_angle

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
