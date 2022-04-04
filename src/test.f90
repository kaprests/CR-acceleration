module tests
    public
contains

subroutine test_spherical_cartesian_coord_change(num_tests_run, num_tests_failed)
    use constants

    implicit none
    integer, intent(inout) :: num_tests_failed, num_tests_run
    double precision :: x, y, z
    double precision :: x2, y2, z2
    double precision :: r, theta, phi
    double precision :: ran0
    integer :: n = 500
    integer :: i

    do i = 1, n, 1
        x = ran0()
        y = ran0()
        z = ran0()
        call cartesian_to_spherical(x, y, z, r, theta, phi)
        call spherical_to_cartesian(r, theta, phi, x2, y2, z2)
        if (abs(x - x2) > 1e-6) then
            print *, "Failed: test_spherical_cartesian_coord_change, line 18"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
        if (abs(y - y2) > 1e-6) then
            print *, "Failed: test_spherical_cartesian_coord_change, line 23"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
        if (abs(z - z2) > 1e-6) then
            print *, "Failed: test_spherical_cartesian_coord_change, line 28"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
    end do

    ! Special case x ==  0 y > 0
    x = 0.d0
    y = 1.d0/sqrt(2.d0)
    z = 1.d0/sqrt(2.d0)
    call cartesian_to_spherical(x, y, z, r, theta, phi)
    call spherical_to_cartesian(r, theta, phi, x2, y2, z2)
    if (abs(x - x2) > 1e-6) then
        print *, "Failed: test_spherical_cartesian_coord_change, line 18"
        num_tests_failed = num_tests_failed + 1
    else if (abs(y - y2) > 1e-6) then
        print *, "Failed: test_spherical_cartesian_coord_change, line 23"
        num_tests_failed = num_tests_failed + 1
    else if (abs(z - z2) > 1e-6) then
        print *, "Failed: test_spherical_cartesian_coord_change, line 28"
        num_tests_failed = num_tests_failed + 1
    end if

    ! Special case x ==  0 y < 0
    x = 0.d0
    y = -1.d0/sqrt(2.d0)
    z = 1.d0/sqrt(2.d0)
    call cartesian_to_spherical(x, y, z, r, theta, phi)
    call spherical_to_cartesian(r, theta, phi, x2, y2, z2)
    if (abs(x - x2) > 1e-6) then
        print *, "Failed: test_spherical_cartesian_coord_change, line 18"
        num_tests_failed = num_tests_failed + 1
    else if (abs(y - y2) > 1e-6) then
        print *, "Failed: test_spherical_cartesian_coord_change, line 23"
        num_tests_failed = num_tests_failed + 1
    else if (abs(z - z2) > 1e-6) then
        print *, "Failed: test_spherical_cartesian_coord_change, line 28"
        num_tests_failed = num_tests_failed + 1
    end if

    num_tests_run = num_tests_run + 1
end subroutine test_spherical_cartesian_coord_change

subroutine test_lorentz_boost(num_tests_run, num_tests_failed)
    implicit none
    integer, intent(inout) :: num_tests_run, num_tests_failed
    double precision :: t, t2, t_prime
    double precision, dimension(3) :: r_vec, r_vec2, r_vec_prime, v_rel_vec
    double precision :: ran0
    integer :: i

    ! initial coords (1, 1, 1, 1)
    t = 1.d0
    r_vec(1) = 1.d0
    r_vec(2) = 1.d0
    r_vec(3) = 1.d0

    ! Boost along x-axis
    v_rel_vec(1) = 0.9d0
    v_rel_vec(2) = 0.0d0
    v_rel_vec(3) = 0.0d0
    call lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    call lorentz_boost(t_prime, r_vec_prime, t2, r_vec2, -v_rel_vec)

    if (abs(t - t2) > 1e-6) then
        print *, "Failed: test_lorentz_boost, x-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(1) - r_vec2(1)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, x-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(2) - r_vec2(2)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, x-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(3) - r_vec2(3)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, x-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if

    ! Boost along y-axis
    v_rel_vec(1) = 0.0d0
    v_rel_vec(2) = 0.9d0
    v_rel_vec(3) = 0.0d0
    call lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    call lorentz_boost(t_prime, r_vec_prime, t2, r_vec2, -v_rel_vec)

    if (abs(t - t2) > 1e-6) then
        print *, "Failed: test_lorentz_boost, y-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(1) - r_vec2(1)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, y-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(2) - r_vec2(2)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, y-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(3) - r_vec2(3)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, y-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if

    ! Boost along z-axis
    t = 1.d0
    v_rel_vec(1) = 0.0d0
    v_rel_vec(2) = 0.0d0
    v_rel_vec(3) = 0.9d0
    t_prime = 0.d0
    r_vec_prime = 0.d0
    t2 = 0.d0
    r_vec2 = 0.d0
    call lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    call lorentz_boost(t_prime, r_vec_prime, t2, r_vec2, -v_rel_vec)

    if (abs(t - t2) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(1) - r_vec2(1)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(2) - r_vec2(2)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(3) - r_vec2(3)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if

    ! Boost along oblique-axis
    t = 1.d0
    r_vec(1) = 1.d0
    r_vec(2) = 1.d0
    r_vec(3) = 1.d0
    v_rel_vec(1) = 0.1d0
    v_rel_vec(2) = 0.2d0
    v_rel_vec(3) = 0.3d0
    t_prime = 0.d0
    r_vec_prime = 0.d0
    t2 = 0.d0
    r_vec2 = 0.d0
    call lorentz_boost0(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    call lorentz_boost0(t_prime, r_vec_prime, t2, r_vec2, -v_rel_vec)

    if (abs(t - t2) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(1) - r_vec2(1)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(2) - r_vec2(2)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    if (abs(r_vec(3) - r_vec2(3)) > 1e-6) then
        print *, "Failed: test_lorentz_boost, z-axis"
        num_tests_failed = num_tests_failed + 1
        return
    end if
    num_tests_run = num_tests_run + 1
end subroutine test_lorentz_boost

subroutine test_rotations(num_tests_run, num_tests_failed)
    use constants

    implicit none
    integer, intent(inout) :: num_tests_run, num_tests_failed
    double precision, dimension(3, 3) :: Ry, Rz, RyRz, RzRy
    double precision, dimension(3) :: z_hat_rot, z_hat_rot_prime, z_hat_rot2
    double precision :: r_unity, theta, phi
    double precision :: ran0
    integer :: i

    do i = 1, 50, 1
        ! z unit vector of rotated in lab frame:
        z_hat_rot = [ran0(), ran0(), ran0()]
        z_hat_rot = z_hat_rot/sqrt(z_hat_rot(1)**2+z_hat_rot(2)**2+z_hat_rot(3)**2)

        ! z unit vector of rotated frame in the rotated frame
        z_hat_rot_prime = [0.d0, 0.d0, 1.d0]

        ! Get rotation angles of prime frame in the lab frame
        call cartesian_to_spherical(z_hat_rot(1), z_hat_rot(2), z_hat_rot(3), r_unity, theta, phi)
        !print *, "z_hat_rot (lab frame): "
        !print *, "x: ", z_hat_rot(1), "y: ", z_hat_rot(2), "z: ", z_hat_rot(3)
        !print *, "r_unity: ", r_unity, "theta: ", theta/pi, "pi ", "phi: ", phi/pi, "pi"

        ! Rotation matrix/matrices
        call euler_Rz(phi, Rz)
        call euler_Ry(theta, Ry)
        RzRy = matmul(Rz, Ry)
        z_hat_rot2 = matmul(Rz, matmul(Ry, z_hat_rot_prime))
        !print *, 'z_hat_rot_prime : ', z_hat_rot_prime
        !print *, 'target: ', z_hat_rot
        !print *, 'back rot: ', matmul(Rz, matmul(Ry, z_hat_rot_prime))
        !print *, 'back rot2: ', matmul(RzRy, z_hat_rot_prime)
        if (abs(z_hat_rot(1) - z_hat_rot2(1)) > 1e-6) then
            print *, "Failed: test_rotations"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
        if (abs(z_hat_rot(2) - z_hat_rot2(2)) > 1e-6) then
            print *, "Failed: test_rotations"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
        if (abs(z_hat_rot(3) - z_hat_rot2(3)) > 1e-6) then
            print *, "Failed: test_rotations"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
    end do

    num_tests_run = num_tests_run + 1
end subroutine test_rotations

subroutine test_small_angle_scattering(num_tests_run, num_tests_failed)
    use constants

    implicit none
    integer, intent(inout) :: num_tests_run, num_tests_failed
    double precision :: theta_max
    double precision :: theta, phi, r_unity
    double precision, dimension(3) :: p_inital, p_scattered
    integer :: i, j

    do i = 1, 10, 1
    theta_max = i * 0.1 * pi
    do j = 1, 500, 1
        ! Isotropically distributed inital angles
        call isotropic_scatter(theta, phi)

        ! Initial direction
        call spherical_to_cartesian(1.d0, theta, phi, p_inital(1), p_inital(2), p_inital(3))

        ! Scatter
        call small_angle_scatter(theta, phi, theta_max)

        ! Scattered direction
        call spherical_to_cartesian(1.d0, theta, phi, p_scattered(1), p_scattered(2), p_scattered(3))
        
        if (acos(dot_product(p_inital, p_scattered)) > theta_max) then
            print *, "Failed: test_small_angle_scattering - Scattering angle > theta_max"
            num_tests_failed = num_tests_failed + 1
            exit
        end if
    end do
    end do

    num_tests_run = num_tests_run + 1
end subroutine test_small_angle_scattering

subroutine test_cubic_spline_small_angle_step_correction(num_tests_run, num_tests_failed)
    use constants, only: pi

    implicit none
    integer, intent(inout) :: num_tests_run, num_tests_failed
    double precision :: cubic_spline_small_angle_step_correction, cs_test
    integer, parameter :: n_theta = 50
    double precision, parameter :: delta_theta = pi/(n_theta)
    integer :: i
    double precision, dimension(n_theta) :: theta_array
    double precision :: stepcorr
    
    do i = 0, n_theta, 1
        theta_array = i * delta_theta
        stepcorr = cubic_spline_small_angle_step_correction(i*delta_theta)
        print *, i*delta_theta, ", ", stepcorr
        !print *, i, " ", i*delta_theta/pi, " ", stepcorr
    end do
    num_tests_run = num_tests_run + 1
end subroutine test_cubic_spline_small_angle_step_correction
    
end module tests

program test
    use constants
    use user_variables
    use tests

    implicit none
    integer :: num_tests_run, num_tests_failed

    num_tests_run = 0
    num_tests_failed = 0
    call test_spherical_cartesian_coord_change(num_tests_run, num_tests_failed)
    call test_lorentz_boost(num_tests_run, num_tests_failed)
    call test_rotations(num_tests_run, num_tests_failed)
    call test_small_angle_scattering(num_tests_run, num_tests_failed)
    call test_cubic_spline_small_angle_step_correction(num_tests_run, num_tests_failed)

    print *, " "
    print *, "TESTS FINISHED"
    print *, "-------------------------------------- "
    print *, "Summary:"
    print *, "Num tests run:      ", num_tests_run
    print *, "Num tests succeded: ", num_tests_run - num_tests_failed
    print *, "Num tests failed:   ", num_tests_failed
    print *, "-------------------------------------- "
end program test
