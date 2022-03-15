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
    use user_variables, only: debug
    use SNR_data, only: t_max; 
    use constants; use particle_data, only: m_p
    use event_internal; use result
    use internal
    use test_var, only: sec, accel

    implicit none
    integer, intent(in) :: set, n_injected
    integer k, n_step
    double precision r, m, f, df, dt, dE, delta, l_0, l_0_0
    double precision r_sh1, r_sh2, phi, theta, phi_v, theta_v, d1, d2, dmax, v_rel
    double precision :: gamma_factor, cos_theta, dt_us
    double precision, dimension(3) :: l_vec, l_vec_us, v_rel_vec
    double precision ran0, R_L, t_shock, v_shock, get_v_rel
    integer, pointer :: pid, A, Z
    double precision, pointer :: E, x(:), t, w
    double precision, dimension(4, 4) :: boost_matrix

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

    do
        ! Stepsize
        df = 1.d-99 ! f_tot_rates(A,Z,E,d1,t)            ! interaction rate (1/yr)
        call scales_charged(m, Z, E, t, w, df, dt, dE)
        l_0 = R_L(E, t)/dble(Z)
        l_0_0 = l_0
        if (l_0 <= 0.d0 .or. dt <= 0.d0) call error('wrong scales', 0)

        ! Step direction
        call isotropic(phi, theta)
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

        ! Convert step from spherical to cartesian coordinates
        ! l_vec: cartesian step in the particle's local region's restframe
        ! l_0: stepsize in particle's local region's frame
        ! theta, phi: particle flight anlges in the local frame
        call spherical_to_cartesian(l_0, theta, phi, l_vec(1), l_vec(2), l_vec(3))

        do k = 1, n_step
            r_sh1 = t_shock(t)                              ! Radial position of shock before step
            d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)          ! Radial position of particle before step

            if (d1 < r_sh1) then
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
                ! TODO: Switch with subroutine call
                call lorentz_boost_matrix(-v_rel_vec, boost_matrix)
                dt_us = dt*boost_matrix(1,1)+l_vec(1)*boost_matrix(1,2)+l_vec(2)*boost_matrix(1,3)+l_vec(3)*boost_matrix(1,4)
                l_vec_us(1) = dt*boost_matrix(2,1)+l_vec(1)*boost_matrix(2,2)+l_vec(2)*boost_matrix(2,3)+l_vec(3)*boost_matrix(2,4)
                l_vec_us(2) = dt*boost_matrix(3,1)+l_vec(1)*boost_matrix(3,2)+l_vec(2)*boost_matrix(3,3)+l_vec(3)*boost_matrix(3,4)
                l_vec_us(3) = dt*boost_matrix(4,1)+l_vec(1)*boost_matrix(4,2)+l_vec(2)*boost_matrix(4,3)+l_vec(3)*boost_matrix(4,4)
                dt = dt_us
                l_vec(1) = l_vec_us(1)
                l_vec(2) = l_vec_us(2)
                l_vec(3) = l_vec_us(3)
            end if
            x(1) = x(1) + l_vec(1)
            x(2) = x(2) + l_vec(2)
            x(3) = x(3) + l_vec(3)

            d2 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)              ! new distance
            t = t + dt
            E = E + dE
            r_sh2 = t_shock(t)

            if (d2 < r_sh2 .and. d1 > r_sh1) then    
                ! we have crossed to the left: US -> DS
                call cartesian_to_spherical(x(1), x(2), x(3), v_rel, theta_v, phi_v)
                v_rel = get_v_rel(v_shock(t))
                if (v_rel <= 0.d0) call error("v_rel <= 0", 0)
                gamma_factor = 1.d0 / (1.d0 - v_rel**2)
                cos_theta = &        
                    cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &    
                    sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &    
                    cos(theta)*cos(theta_v)
                E = gamma_factor * E * (1.d0 - v_rel*cos_theta)
                accel = 1
            else if (d2 > r_sh2 .and. d1 < r_sh1) then         
                ! we have crossed to the right: DS -> US
                call cartesian_to_spherical(x(1), x(2), x(3), v_rel, theta_v, phi_v)
                v_rel = get_v_rel(v_shock(t))
                if (v_rel <= 0.d0) call error("v_rel <= 0", 0)
                gamma_factor = 1.d0 / (1.d0 - v_rel**2)
                cos_theta = &        
                    cos(phi)*sin(theta)*cos(phi_v)*sin(theta_v) + &    
                    sin(phi)*sin(theta)*sin(phi_v)*sin(theta_v) + &    
                    cos(theta)*cos(theta_v)
                E = gamma_factor * E * (1.d0 + v_rel*cos_theta)
                accel = 1
            end if

            v_rel = 0.75d0*v_shock(t)
            dmax = 3.d0*l_0_0/v_rel
            f = f + df*dt                      ! \int dt f(t)
            delta = exp(-f)                    ! exp(-\int dt f(t))

            ! exit, if a) too late, b) too far down-stream, or c) scattering:
            if (t > t_max .or. d2 < r_sh2 - dmax .or. r > delta) then
                if (t > t_max .or. d2 < r_sh2 - dmax) then  ! we're tired or trapped behind
                    !              write(*,*) 'tired',n_in,n_out
                    call store(pid, E, w)
                    print *, "Exit energy: ", E
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
    write(*,*) 'store: ',pid,i
end subroutine store


end module random_walk
