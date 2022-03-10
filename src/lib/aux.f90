! File: aux.f90
! Auxillary/utility subroutines and functions


subroutine error(string, s)
    ! error handling
    character(len=*), intent(in) :: string
    integer, intent(in) :: s
    integer, save :: n_warn
    integer :: n_warn_max = 100

    if (s == 1 .or. s == 11) then                   ! warning message
        write (*, *)
        write (*, *) 'Warning:'
        write (*, *) string
        write (99, *) string
        if (s == 1) then                         ! warning
            n_warn = n_warn + 1
            if (n_warn > n_warn_max) then
                write (*, *)
                write (*, *) 'more than', n_warn_max, ' warnings!'
                write (*, *)
                write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write (*, *) '!  ELMAG 3.01 stops program excecution  !'
                write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                stop
            end if
        end if
    end if
    if (s == 0) then                    ! error
        write (*, *)
        write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write (*, *) '!   a serious error:                    !'
        write (99, *) '!   a serious error:                    !'
        write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write (*, *)
        write (*, *) string
        write (99, *) string
        write (*, *)
        write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write (*, *) '!  ELMAG 3.01 stops program excecution  !'
        write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        stop
    end if
end subroutine error


function ran0()
    ! random number generator from Numerical Recipes (Fortran90)
    use internal, only: iseed

    implicit none
    integer, parameter :: K4B = selected_int_kind(9)
    !integer(K4B), intent(inout) :: iseed
    double precision ran0
    integer(K4B), parameter :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
    real, save :: am
    integer(K4B), save :: ix = -1, iy = -1, k

    if (iseed <= 0 .or. iy < 0) then
        am = nearest(1.e0, -1.e0)/IM
        iy = ior(ieor(888889999, abs(iseed)), 1)
        ix = ieor(777755555, abs(iseed))
        iseed = abs(iseed) + 1
    end if
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy/IQ
    iy = IA*(iy - k*IQ) - IR*k
    if (iy < 0) iy = iy + IM
    ran0 = am*ior(iand(IM, ieor(ix, iy)), 1)
end function ran0


subroutine set_theta_max(theta)
    ! Sets default/user provded theta max
    use user_variables, only: theta_max, theta_max_set
    implicit none
    double precision, intent(inout) :: theta
    theta = theta_max
end subroutine set_theta_max


subroutine isotropic(theta, phi)
    use constants

    implicit none
    double precision phi, theta, r, x
    double precision ran0

    r = ran0()
    x = -1.d0 + 2.d0*r                      ! generate x=cos(theta)
    theta = acos(x)
    r = ran0()
    phi = two_pi*r
end subroutine isotropic


subroutine scattering_angle(theta, phi, theta_max)
    ! Random small angle within a cone centered around z-axis
    use constants, only: pi, two_pi

    implicit none
    double precision, intent(inout) :: phi, theta
    double precision, intent(in) :: theta_max
    double precision :: ran0, z

    ! Random angle within the max scattering cone
    z = cos(theta_max) + (1 - cos(theta_max))*ran0()
    theta = acos(z) ! Theta within max
    phi = two_pi*ran0() ! Azimuthal angle phi isotropic
end subroutine scattering_angle


! Not in use
subroutine max_scattering_angle(theta_max_computed, v_shock, E_particle)
    ! Computes the loss cone angle, and sets max_pitch scattering angle
    ! to some fraction of cone angle.
    ! Currently not is use
    use constants, only: pi
    use particle_data, only: m_p
    use user_variables, only: theta_max, theta_max_set

    implicit none
    double precision, intent(out) :: theta_max_computed
    double precision, intent(in) :: v_shock, E_particle
    double precision :: v_particle, v_p
    double precision :: cos_theta_cone, theta_cone

    if (theta_max_set) then
        ! Use user provided theta max
        theta_max_computed = theta_max
    else
        v_p = v_particle(E_particle, m_p)
        if (v_shock > v_p) then
            ! isotropic -- E.g. particles injected in front of UR shock (before overtaken 1st time)
            !theta_max_computed = pi
            theta_max_computed = 1/(1/sqrt(1 - v_shock**2)) !approx for UR shocks from Achterberg et at.
        else
            ! Compute loss cone opening, theta_cone
            ! Set max scattering, theta_max, to 100% of loss cone angle
            cos_theta_cone = v_shock/v_p
            if (abs(cos_theta_cone) > 1) then
                print *, "v_shock: ", v_shock
                print *, "E: ", E_particle
                print *, "m: ", m_p
                print *, "cos_theta_cone: ", cos_theta_cone
                call error("cosine exceeds 1, max_scattering_angle", 0)
            end if
            theta_cone = acos(cos_theta_cone)
            theta_max_computed = 1.0*theta_cone
        end if
    end if
end subroutine max_scattering_angle


subroutine euler_RyRz(theta, phi, R)
    implicit none
    double precision, intent(in) :: theta, phi
    double precision, intent(inout) :: R(3, 3)
    double precision :: ct, cp, st, sp

    ct = cos(theta)
    cp = cos(phi)
    st = sin(theta)
    sp = sin(phi)

    ! R(column, row)
    R(1, 1) = ct*cp
    R(1, 2) = sp
    R(1, 3) = -st*cp

    R(2, 1) = -ct*sp
    R(2, 2) = cp
    R(2, 3) = st*sp

    R(3, 1) = st
    R(3, 2) = 0.d0
    R(3, 3) = ct
end subroutine euler_RyRz


subroutine radially_outward(theta_rad, phi_rad, x1, x2, x3)
    ! Finds the angles corresponding to radially out at point x, i.e. shock normal
    use constants, only: pi, two_pi

    implicit none
    double precision, intent(in) :: x1, x2, x3
    double precision, intent(inout) :: phi_rad, theta_rad
    theta_rad = atan2(sqrt(x1**2 + x2**2), x3)
    phi_rad = atan(x2/x1)
    if (x1 == 0) then
        if (x2 > 0) then
            phi_rad = pi/4
        else if (x2 < 0) then
            phi_rad = 3*pi/4
        else
            phi_rad = pi
        end if
    end if
    if (x1 < 0.d0 .and. x2 > 0) phi_rad = phi_rad + pi
    if (x1 < 0.d0 .and. x2 < 0) phi_rad = phi_rad + pi
    if (x1 > 0.d0 .and. x2 < 0) phi_rad = phi_rad + two_pi
end subroutine radially_outward


double precision function get_v_2(v_shock) result(v)
    ! Move to functions.f90
    use user_variables, only: gamma_sh

    implicit none
    double precision, intent(in) :: v_shock
    if (v_shock <= 0.1 .and. 0 < v_shock) then
        ! Non relativistic regime
        v = 0.75d0*v_shock
    else if (0.9 <= v_shock .and. v_shock < 1) then
        ! (Ultra) relativistic regime
        v = 1.d0 - 1.d0/(gamma_sh**2)
    else if (0.1 < v_shock .and. v_shock < 0.9) then
        ! In between regimes
        v = 0.75d0*v_shock
        call error("v_shock in between regimes 0.1 < v_shock < 0.9", 1)
    else
        call error("v_shock outside of allowed range [0, 1]", 0)
    end if
end function get_v_2


double precision function spacetime_interval(t, r_vec)
    ! "norm" of four vector
    ! metric = diag(-1, 1, 1, 1)
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(3), intent(in) :: r_vec

    spacetime_interval = -1*t**2 + dot_product(r_vec, r_vec)
end function spacetime_interval


subroutine spherical_to_cartesian(r, theta, phi, x, y, z)
    ! Convert from spherical to cartesian coordinates
    implicit none
    double precision, intent(in) :: r, phi, theta
    double precision, intent(out) :: x, y, z

    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
end subroutine spherical_to_cartesian


subroutine parallel_projection(v, u, v_parallel)
    implicit none
    double precision, dimension(3), intent(in) :: v, u
    double precision, dimension(3), intent(out) :: v_parallel

    v_parallel = (dot_product(u, v) / dot_product(u, u)) * u
end subroutine parallel_projection


subroutine orthogonal_projection(v, u, v_ortogonal)
    implicit none
    double precision, dimension(3), intent(in) :: v, u
    double precision, dimension(3), intent(out) :: v_ortogonal
    double precision, dimension(3) :: v_parallel

    call parallel_projection(v, u, v_parallel)
    v_ortogonal = v - v_parallel
end subroutine orthogonal_projection


double precision function v_prime(v, v_rel)
    ! Velocity boost transformation
    implicit none
    double precision, intent(in) :: v, v_rel

    v_prime = (v - v_rel)/(1 - v*v_rel)
end function v_prime


subroutine lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    ! Lorentz boost
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(3), intent(in) :: r_vec, v_rel_vec
    double precision, intent(out) :: t_prime
    double precision, dimension(3), intent(out) :: r_vec_prime
    double precision :: gamma_factor
    double precision, dimension(3) :: r_parallel, r_orthogonal, r_parallel_prime
    double precision, dimension(4) :: four_vec_parallel, four_vec_parallel_prime
    double precision, dimension(4, 4) :: boost_matrix
    integer :: i, j
    double precision :: spacetime_interval

    gamma_factor = 1.0/sqrt(1.0 - dot_product(v_rel_vec, v_rel_vec))
    call parallel_projection(r_vec, v_rel_vec, r_parallel)
    call orthogonal_projection(r_vec, v_rel_vec, r_orthogonal)
    four_vec_parallel = [t, r_parallel(1), r_parallel(2), r_parallel(3)]
    boost_matrix = transpose(reshape([&
        gamma_factor, -gamma_factor*v_rel_vec(1), -gamma_factor*v_rel_vec(2), -gamma_factor*v_rel_vec(3), &
        -gamma_factor*v_rel_vec(1), gamma_factor, 0.d0, 0.d0, &
        -gamma_factor*v_rel_vec(2), 0.d0, gamma_factor, 0.d0, &
        -gamma_factor*v_rel_vec(3), 0.d0, 0.d0, gamma_factor  &
        ], shape(boost_matrix)&
    ))
    four_vec_parallel_prime = 0.d0
    do i = 1, 4, 1 ! column i
    do j = 1, 4, 1 ! row j
        four_vec_parallel_prime(i) = four_vec_parallel_prime(i) + boost_matrix(j, i) * four_vec_parallel(i)
    end do
    end do
    t_prime = four_vec_parallel_prime(1)
    r_parallel_prime = four_vec_parallel_prime(2:)
    r_vec_prime = r_parallel_prime + r_orthogonal
    if (spacetime_interval(t, r_vec) - spacetime_interval(t_prime, r_vec_prime) > 0.0001) then
        print *, "-------------------------"
        print *, "t: ", t
        print *, "r_vec: ", r_vec
        print *, "t_prime: ", t_prime
        print *, "r_vec_prime: ", r_vec_prime
        print *, "Delta ds^2: ", spacetime_interval(t, r_vec) - spacetime_interval(t_prime, r_vec_prime)
        call error("Erroneous boost, spacetime interval not invariant.", 0)
    end if
end subroutine lorentz_boost


subroutine upstream_to_downstream_boost(t_us, r_us, t_ds, r_ds, v_shock_us, pos_vec_us)
    ! Upstream restframe to downstream restframe
    implicit none
    double precision, intent(in) :: t_us, v_shock_us
    double precision, dimension(3), intent(in) :: r_us, pos_vec_us
    double precision, intent(out) :: t_ds
    double precision, dimension(3), intent(out) :: r_ds
    double precision :: v_rel, theta, phi
    double precision, dimension(3) :: v_rel_vec
    double precision :: get_v_2

    v_rel = get_v_2(v_shock_us)
    call radially_outward(theta, phi, pos_vec_us(1), pos_vec_us(2), pos_vec_us(3))
    call spherical_to_cartesian(v_rel, theta, phi, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
    call lorentz_boost(t_us, r_us, t_ds, r_ds, v_rel_vec)
end subroutine upstream_to_downstream_boost


subroutine downstream_to_upstream_boost(t_ds, r_ds, t_us, r_us, v_shock_us, pos_vec_us)
    ! Downstream restframe to upstream restframe boost
    implicit none
    double precision, intent(in) :: t_ds, v_shock_us
    double precision, dimension(3), intent(in) :: r_ds, pos_vec_us
    double precision, intent(out) :: t_us
    double precision, dimension(3), intent(out) :: r_us
    double precision :: v_rel, theta, phi
    double precision, dimension(3) :: v_rel_vec
    double precision :: get_v_2

    v_rel = get_v_2(v_shock_us)
    call radially_outward(theta, phi, pos_vec_us(1), pos_vec_us(2), pos_vec_us(3))
    call spherical_to_cartesian(v_rel, theta, phi, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
    v_rel_vec = -1.d0 * v_rel_vec ! Minus sign since US flows radially inward
    call lorentz_boost(t_ds, r_ds, t_us, r_us, v_rel_vec)
end subroutine downstream_to_upstream_boost


subroutine upstream_to_shockfront_boost(t_us, r_us, t_sh, r_sh, v_shock_us, pos_vec_us)
    ! Upstream to shockfront boost
    implicit none
    double precision, intent(in) :: t_us, v_shock_us
    double precision, dimension(3), intent(in) :: r_us, pos_vec_us
    double precision, intent(out) :: t_sh
    double precision, dimension(3), intent(out) :: r_sh
    double precision :: theta, phi
    double precision, dimension(3) :: v_shock_vec

    call radially_outward(theta, phi, pos_vec_us(1), pos_vec_us(2), pos_vec_us(3))
    call spherical_to_cartesian(v_shock_us, theta, phi, v_shock_vec(1), v_shock_vec(2), v_shock_vec(3))
    call lorentz_boost(t_us, r_us, t_sh, r_sh, v_shock_vec)
end subroutine upstream_to_shockfront_boost


!subroutine lorentz_boost_time_component(t, t_prime, v_rel)
!    implicit none
!    double precision, intent(in) :: t, v_rel
!    double precision, intent(out) :: t_prime
!    double precision :: gamma_factor
!    gamma_factor = 1.d0/sqrt(1.d0 - v_rel**2)
!    t_prime = gamma_factor * ()
!end subroutine lorentz_boost_time_component









