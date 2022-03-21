! File: aux.f90
! Auxillary/utility subroutines and functions

subroutine error(string, s)
    !Error handling
    implicit none
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
    !random number generator from Numerical Recipes (Fortran90)
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

subroutine isotropic(phi, theta)
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

subroutine euler_Ry(theta, R)
    implicit none
    double precision, intent(in) :: theta
    double precision, dimension(3,3), intent(out) :: R
    double precision :: ct, st

    ct = cos(theta)
    st = sin(theta)

    ! First column
    R(1, 1) = ct
    R(2, 1) = 0
    R(3, 1) = -st

    ! Second column
    R(1, 2) = 0
    R(2, 2) = 1
    R(3, 2) = 0
    
    ! Third column
    R(1, 3) = st
    R(2, 3) = 0
    R(3, 3) = ct
end subroutine euler_Ry

subroutine euler_Rz(theta, R)
    implicit none
    double precision, intent(in) :: theta
    double precision, dimension(3,3), intent(out) :: R
    double precision :: ct, st

    ct = cos(theta)
    st = sin(theta)

    ! First column
    R(1, 1) = ct
    R(2, 1) = st
    R(3, 1) = 0

    ! Second column
    R(1, 2) = -st
    R(2, 2) = ct
    R(3, 2) = 0
    
    ! Third column
    R(1, 3) = 0
    R(2, 3) = 0
    R(3, 3) = 1
end subroutine euler_Rz

subroutine euler_RyRz(theta, phi, RyRz)
    implicit none
    double precision, intent(in) :: theta, phi
    double precision, dimension(3, 3), intent(inout) :: RyRz
    double precision, dimension(3, 3) :: Ry, Rz

    call euler_Ry(theta, Ry)
    call euler_Rz(theta, Rz)
    RyRz = matmul(Ry, Rz)
end subroutine euler_RyRz

subroutine euler_RyRz2(theta, phi, R)
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
    R(2, 1) = sp
    R(3, 1) = -st*cp

    R(1, 2) = -ct*sp
    R(2, 2) = cp
    R(3, 2) = st*sp

    R(1, 3) = st
    R(2, 3) = 0.d0
    R(3, 3) = ct
end subroutine euler_RyRz2

double precision function get_v_rel(v_shock)
    implicit none
    double precision, intent(in) :: v_shock
    double precision :: gamma_sh

    if (v_shock <= 0.1 .and. 0 < v_shock) then
        ! Non relativistic regime
        get_v_rel = 0.75d0*v_shock
    else if (0.9 <= v_shock .and. v_shock < 1) then
        ! (Ultra) relativistic regime
        gamma_sh = 1.d0/sqrt(1.d0 - v_shock**2)
        get_v_rel = 1.d0 - 1.d0/(gamma_sh**2)
    else if (0.1 < v_shock .and. v_shock < 0.9) then
        ! In between regimes
        get_v_rel = 0.75d0*v_shock
        call error("v_shock in between regimes 0.1 < v_shock < 0.9", 1)
    else
        call error("v_shock outside of allowed range [0, 1]", 0)
    end if
end function get_v_rel

subroutine cartesian_to_spherical(x, y, z, r, theta, phi)
    use constants, only: pi, two_pi

    implicit none
    double precision, intent(in) :: x, y, z
    double precision, intent(out) :: r, theta, phi

    r = sqrt(x**2 + y**2 + z**2)
    theta = atan2(sqrt(x**2 + y**2), z)
    phi = atan(y/x)
    if (x == 0) then
        if (y > 0) then
            phi = pi/4
        else if (y < 0) then
            phi = 3*pi/4
        else
            phi = pi
        end if
    end if
    if (x < 0.d0 .and. y > 0) phi = phi + pi
    if (x < 0.d0 .and. y < 0) phi = phi + pi
    if (x > 0.d0 .and. y < 0) phi = phi + two_pi
end subroutine cartesian_to_spherical

subroutine spherical_to_cartesian(r, theta, phi, x, y, z)
    implicit none
    double precision, intent(in) :: r, theta, phi
    double precision, intent(out) :: x, y, z

    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)
end subroutine spherical_to_cartesian

subroutine parallel_projection(v, u, v_parallel)
    ! parallel vector projection
    implicit none
    double precision, dimension(3), intent(in) :: v, u
    double precision, dimension(3), intent(out) :: v_parallel

    v_parallel = (dot_product(u, v)/dot_product(u, u))*u
end subroutine parallel_projection

subroutine orthogonal_projection(v, u, v_orthogonal)
    implicit none
    double precision, dimension(3), intent(in) :: v, u
    double precision, dimension(3), intent(out) :: v_orthogonal
    double precision, dimension(3) :: v_parallel

    call parallel_projection(v, u, v_parallel)
    v_orthogonal = v - v_parallel
end subroutine orthogonal_projection

double precision function v_prime(v, v_rel)
    ! Velocity boost tranform
    implicit none
    double precision, intent(in) :: v, v_rel

    v_prime = (v - v_rel)/(1 - v*v_rel)
end function v_prime

double precision function spacetime_interval(t, r_vec)
    ! "norm" of four vector
    ! metric = diag(-1, 1, 1, 1)
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(3), intent(in) :: r_vec

    spacetime_interval = -1*t**2 + dot_product(r_vec, r_vec)
end function spacetime_interval

subroutine lorentz_boost0(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    ! Lorentz boost
    ! Parameters:
    !   t: time component in original frame
    !   r_vec: spatial component in original frame
    !   t_prime: time component in the boosted frame
    !   r_vec_prime: spatial component on boosted frame
    !   v_rel_vec: relative velocity vector of the two frames
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

    gamma_factor = 1.d0/sqrt(1.d0 - dot_product(v_rel_vec, v_rel_vec))
    call parallel_projection(r_vec, v_rel_vec, r_parallel)
    call orthogonal_projection(r_vec, v_rel_vec, r_orthogonal)
    four_vec_parallel = [t, r_parallel(1), r_parallel(2), r_parallel(3)]
    boost_matrix = transpose(reshape([ &
        ! Row1
        gamma_factor, -gamma_factor*v_rel_vec(1), &
        -gamma_factor*v_rel_vec(2), -gamma_factor*v_rel_vec(3), &
        ! Row2
        -gamma_factor*v_rel_vec(1), gamma_factor, 0.d0, 0.d0, &
        ! Row3
        -gamma_factor*v_rel_vec(2), 0.d0, gamma_factor, 0.d0, &
        ! Row4
        -gamma_factor*v_rel_vec(3), 0.d0, 0.d0, gamma_factor &
        ], shape(boost_matrix) &
        ))
    four_vec_parallel_prime = matmul(boost_matrix, four_vec_parallel)
    t_prime = four_vec_parallel_prime(1)
    r_parallel_prime = four_vec_parallel_prime(2:)
    r_vec_prime = r_parallel_prime + r_orthogonal
    if (spacetime_interval(t, r_vec) - spacetime_interval(t_prime, r_vec_prime) > 0.0001) then
        print *, "-------------------------"
        print *, "t: ", t
        print *, "r_vec: ", r_vec
        print *, "t_prime: ", t_prime
        print *, "r_vec_prime: ", r_vec_prime
        print *, "Delta ds^2: ", spacetime_interval(t, r_vec)-spacetime_interval(t_prime, r_vec_prime)
        call error("Erroneous boost, spacetime interval not invariant.", 0)
    end if
end subroutine lorentz_boost0

subroutine lorentz_boost_matrix(v_rel_vec, boost_matrix)
    implicit none
    double precision, dimension(3), intent(in) :: v_rel_vec
    double precision, dimension(4, 4), intent(out) :: boost_matrix
    double precision :: gamma_factor, v_squared, vx, vy, vz

    vx = v_rel_vec(1)
    vy = v_rel_vec(2)
    vz = v_rel_vec(3)
    v_squared = dot_product(v_rel_vec, v_rel_vec)
    gamma_factor = 1.d0/sqrt(1.d0 - v_squared)
    boost_matrix = 0.d0
    boost_matrix = transpose(reshape([ &
        ! row1
        gamma_factor, &                         ! (1, 1)
        -gamma_factor*vx, &                   ! (1, 2)
        -gamma_factor*vy, &                   ! (1, 3)
        -gamma_factor*vz, &                   ! (1, 4)
        ! row2
        -gamma_factor*vx, &                                         ! (2, 1)
        1.d0 + (gamma_factor - 1.d0)*((vx*vx)/(v_squared)), &       ! (2, 2)
        (gamma_factor - 1.d0)*((vx*vy)/(v_squared)), &              ! (2, 3)
        (gamma_factor - 1.d0)*((vx*vz)/(v_squared)), &              ! (2, 4)
        ! row3
        -gamma_factor*vy, &                                         ! (3, 1)
        (gamma_factor - 1.d0)*((vy*vx)/(v_squared)), &              ! (3, 2)
        1.d0 + (gamma_factor - 1.d0)*((vy*vy)/(v_squared)), &       ! (3, 3)
        (gamma_factor - 1.d0)*((vy*vz)/(v_squared)), &              ! (3, 4)
        ! row4
        -gamma_factor*vz, &                                         ! (4, 1)
        (gamma_factor - 1.d0)*((vz*vx)/(v_squared)), &              ! (4, 2)
        (gamma_factor - 1.d0)*((vz*vy)/(v_squared)), &              ! (3, 2)
        1.d0 + (gamma_factor - 1.d0)*((vz*vz)/(v_squared)) &        ! (4, 4)
        ], shape(boost_matrix)))
    !print *, "---------------------------------"
    !print *, v_rel_vec
    !print *, boost_matrix(1, :)
    !print *, boost_matrix(2, :)
    !print *, boost_matrix(3, :)
    !print *, boost_matrix(4, :)
    !print *, "---------------------------------"
end subroutine lorentz_boost_matrix

subroutine lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(3), intent(in) :: r_vec, v_rel_vec
    double precision, intent(out) :: t_prime
    double precision, dimension(3), intent(out) :: r_vec_prime
    double precision :: gamma_factor
    double precision, dimension(4) :: four_vec, four_vec_prime
    double precision, dimension(4, 4) :: boost_matrix
    integer :: i, j
    double precision :: spacetime_interval
    double precision :: vrel_squared

    call lorentz_boost_matrix(v_rel_vec, boost_matrix)
    four_vec = [t, r_vec(1), r_vec(2), r_vec(3)]
    four_vec_prime = matmul(boost_matrix, four_vec)
    !do i = 1, 4, 1
    !do j = 1, 4, 1
    !    four_vec_prime(i) = four_vec_prime(i) + boost_matrix(j, i)*four_vec(i)
    !end do
    !end do
    !four_vec_prime = &
    !    four_vec(1) * boost_matrix(1, :) + &
    !    four_vec(2) * boost_matrix(2, :) + &
    !    four_vec(3) * boost_matrix(3, :) + &
    !    four_vec(4) * boost_matrix(4, :)
    t_prime = four_vec_prime(1)
    r_vec_prime = four_vec_prime(2:)
end subroutine lorentz_boost

subroutine upstream_to_downstream_boost(t_us, r_us, t_ds, r_ds, v_shock_us, pos_vec_us)
    ! Upstream restframe to downstream restframe
    implicit none
    double precision, intent(in) :: t_us, v_shock_us
    double precision, dimension(3), intent(in) :: r_us, pos_vec_us
    double precision, intent(out) :: t_ds
    double precision, dimension(3), intent(out) :: r_ds
    double precision :: v_rel, r, theta, phi
    double precision, dimension(3) :: v_rel_vec
    double precision :: get_v_rel

    v_rel = get_v_rel(v_shock_us)
    call cartesian_to_spherical(pos_vec_us(1), pos_vec_us(2), pos_vec_us(3), r, theta, phi)
    call spherical_to_cartesian(v_rel, theta, phi, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
    call lorentz_boost(t_us, r_us, t_ds, r_ds, v_rel_vec)
end subroutine upstream_to_downstream_boost

subroutine downstream_to_upstream_boost(t_ds, r_ds, t_us, r_us, v_shock_us, time_us, pos_vec_us)
    ! Downstream restframe to upstream restframe boost
    implicit none
    double precision, intent(in) :: t_ds, v_shock_us, time_us
    double precision, dimension(3), intent(in) :: r_ds, pos_vec_us
    double precision, intent(out) :: t_us
    double precision, dimension(3), intent(out) :: r_us
    double precision :: v_rel, r, theta, phi
    double precision, dimension(3) :: v_rel_vec
    double precision :: get_v_rel
    double precision :: time_ds
    double precision, dimension(3) :: pos_vec_ds

    ! Transform pos_vec_us to pos_vec_ds
    call upstream_to_downstream_boost(time_us, pos_vec_us, time_ds, pos_vec_ds, v_shock_us, pos_vec_us)
    v_rel = get_v_rel(v_shock_us)
    call cartesian_to_spherical(pos_vec_ds(1), pos_vec_ds(2), pos_vec_ds(3), r, theta, phi)
    call spherical_to_cartesian(v_rel, theta, phi, v_rel_vec(1), v_rel_vec(2), v_rel_vec(3))
    v_rel_vec = -1.d0*v_rel_vec ! Minus sign since US flows radially inward
    call lorentz_boost(t_ds, r_ds, t_us, r_us, v_rel_vec)
end subroutine downstream_to_upstream_boost

subroutine upstream_to_shockfront_boost(t_us, r_us, t_sh, r_sh, v_shock_us, pos_vec_us)
    ! Upstream to shockfront boost
    implicit none
    double precision, intent(in) :: t_us, v_shock_us
    double precision, dimension(3), intent(in) :: r_us, pos_vec_us
    double precision, intent(out) :: t_sh
    double precision, dimension(3), intent(out) :: r_sh
    double precision :: r, theta, phi
    double precision, dimension(3) :: v_shock_vec

    call cartesian_to_spherical(pos_vec_us(1), pos_vec_us(2), pos_vec_us(3), r, theta, phi)
    call spherical_to_cartesian(v_shock_us, theta, phi, v_shock_vec(1), v_shock_vec(2), v_shock_vec(3))
    call lorentz_boost(t_us, r_us, t_sh, r_sh, v_shock_vec)
end subroutine upstream_to_shockfront_boost
