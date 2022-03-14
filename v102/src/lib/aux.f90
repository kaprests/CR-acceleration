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


subroutine cartesian_to_spherical_vec(x_vec, r_vec)
    ! Converts from cartesian to spheical coordinates
    ! The angles are then defining the radially outward direction
    ! at the given position
    use constants, only: pi, two_pi

    implicit none
    double precision, dimension(3), intent(in) :: x_vec
    double precision, dimension(3), intent(out) :: r_vec
    double precision :: r, theta, phi

    r = sqrt(x_vec(1)**2 + x_vec(2)**2 + x_vec(3)**2)
    theta = atan2(sqrt(x_vec(1)**2 + x_vec(2)**2), x_vec(3))
    phi = atan(x_vec(2)/x_vec(1))
    if (x_vec(1) == 0) then
        if (x_vec(2) > 0) then
            phi = pi/4
        else if (x_vec(2) < 0) then
            phi = pi
        end if
    end if
    if (x_vec(1) < 0.d0 .and. x_vec(2) > 0) phi = phi + pi
    if (x_vec(1) < 0.d0 .and. x_vec(2) < 0) phi = phi + pi
    if (x_vec(1) > 0.d0 .and. x_vec(2) < 0) phi = phi + two_pi
    r(1) = r
    r(2) = theta
    r(3) = phi
end subroutine cartesian_to_spherical_vec


subroutine cartesian_pos_to_spherical_angles(x_vec, theta, phi)
    implicit none
    double precision, dimension(3), intent(in) :: x_vec
    double precision, intent(out) :: theta, phi
    double precision, dimension(3) :: r_vec
    
    call cartesian_to_spherical(x_vec, r_vec)
    theta = r_vec(2)
    phi = r_vec(3)
end subroutine cartesian_pos_to_spherical_angles


subroutine spherical_to_cartesian_vec(r_vec, x_vec)
    ! Converts from spherical to cartesian coordinates
    implicit none
    double precision, dimension(3), intent(in) :: r_vec
    double precision, dimension(3), intent(out) :: x_vec

    x(1) = r_vec(1) * sin(r_vec(2)) * cos(r_vec(3))
    x(2) = r_vec(1) * sin(r_vec(2)) * sin(r_vec(3))
    x(3) = r_vec(1) * cos(r_vec(2))
end subroutine spherical_to_cartesian_vec


subroutine spherical_to_cartesian(r, theta, phi, x, y, z)
    implicit none
    double precision, dimension(3), intent(in) :: r_vec
    double precision, intent(out) :: x, y, z
end subroutine spherical_to_cartesian


subroutine parallel_projection(v, u, v_parallel)
    ! parallel vector projection
    implicit none
    double precision, dimension(3), intent(in) :: v, u
    double precision, dimension(3), intent(out) :: v_parallel

    v_parallel = (dot_product(u, v) / dot_product(u, u)) * u
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


subroutine lortentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
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
    
    gamma_factor = 1.0/sqrt(1.0 - dot_product(v_rel_vec, v_rel_vec))    
    call parallel_projection(r_vec, v_rel_vec, r_parallel)    
    call orthogonal_projection(r_vec, v_rel_vec, r_orthogonal)    
    four_vec_parallel = [t, r_parallel(1), r_parallel(2), r_parallel(3)]    
    boost_matrix = transpose(reshape([&    
        ! Row1
        gamma_factor, -gamma_factor*v_rel_vec(1), &
        -gamma_factor*v_rel_vec(2), -gamma_factor*v_rel_vec(3), &
        ! Row2
        -gamma_factor*v_rel_vec(1), gamma_factor, 0.d0, 0.d0, &    
        ! Row3
        -gamma_factor*v_rel_vec(2), 0.d0, gamma_factor, 0.d0, &    
        ! Row4
        -gamma_factor*v_rel_vec(3), 0.d0, 0.d0, gamma_factor  &    
        ], shape(boost_matrix)&    
    ))    
    four_vec_parallel_prime = 0.d0    
    do i = 1, 4, 1 ! column i    
    do j = 1, 4, 1 ! row j    
        four_vec_parallel_prime(i) = &
            four_vec_parallel_prime(i) + boost_matrix(j, i) * four_vec_parallel(i)
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
        print *,"Delta ds^2: ",spacetime_interval(t, r_vec)-spacetime_interval(t_prime, r_vec_prime)
        call error("Erroneous boost, spacetime interval not invariant.", 0)
    end if
end subroutine lortentz_boost






















