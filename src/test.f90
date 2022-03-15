subroutine test_spherical_cartesian_coord_change(num_tests_run, num_tests_failed)
    use constants

    implicit none
    integer, intent(inout) :: num_tests_failed, num_tests_run
    double precision :: x, y, z
    double precision :: x2, y2, z2
    double precision :: r, theta, phi
    double precision :: ran0
    integer :: n = 100
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

    print *, "BOOST x-ax: "
    print *, "---------------------------"
    print *, "v_rel: ", v_rel_vec
    print *, "Original frame: "
    print *, "t: ", t
    print *, "r_vec: ", r_vec
    print *, ""
    print *, "Boosted frame: "
    print *, "t_prime: ", t_prime
    print *, "r_vec_prime: ", r_vec_prime
    print *, ""
    print *, "t2: ", t2
    print *, "r_vec2: ", r_vec2
    print *, "---------------------------"
    print *, ""

    ! Boost along y-axis
    v_rel_vec(1) = 0.0d0
    v_rel_vec(2) = 0.9d0
    v_rel_vec(3) = 0.0d0
    call lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    call lorentz_boost(t_prime, r_vec_prime, t2, r_vec2, -v_rel_vec)

    print *, "BOOST y-ax: "
    print *, "---------------------------"
    print *, "v_rel: ", v_rel_vec
    print *, "Original frame: "
    print *, "t: ", t
    print *, "r_vec: ", r_vec
    print *, ""
    print *, "Boosted frame: "
    print *, "t_prime: ", t_prime
    print *, "r_vec_prime: ", r_vec_prime
    print *, ""
    print *, "t2: ", t2
    print *, "r_vec2: ", r_vec2
    print *, "---------------------------"
    print *, ""

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
    call lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    call lorentz_boost(t_prime, r_vec_prime, t2, r_vec2, -v_rel_vec)

    print *, "BOOST oblique-ax: "
    print *, "---------------------------"
    print *, "v_rel: ", v_rel_vec
    print *, "Original frame: "
    print *, "t: ", t
    print *, "r_vec: ", r_vec
    print *, ""
    print *, "Boosted frame: "
    print *, "t_prime: ", t_prime
    print *, "r_vec_prime: ", r_vec_prime
    print *, ""
    print *, "t2: ", t2
    print *, "r_vec2: ", r_vec2
    print *, "---------------------------"
    print *, ""

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

    print *, "BOOST oblique-ax: "
    print *, "---------------------------"
    print *, "v_rel: ", v_rel_vec
    print *, "Original frame: "
    print *, "t: ", t
    print *, "r_vec: ", r_vec
    print *, ""
    print *, "Boosted frame: "
    print *, "t_prime: ", t_prime
    print *, "r_vec_prime: ", r_vec_prime
    print *, ""
    print *, "t2: ", t2
    print *, "r_vec2: ", r_vec2
    print *, "---------------------------"
    print *, ""

    ! Low boost velocity check (vs gallilean)
    print *, "LOW VELOCITY CHECK: "
    v_rel_vec(1) = 0.01d0
    v_rel_vec(2) = 0.0d0
    v_rel_vec(3) = 0.0d0
    t_prime = 0.d0
    r_vec_prime = 0.d0
    t2 = 0.d0
    r_vec2 = 0.d0
    call lorentz_boost(t, r_vec, t_prime, r_vec_prime, v_rel_vec)
    print *, "Lorentz boosted: "
    print *, "t_prime: ", t_prime
    print *, "r_vec_prime: ", r_vec_prime
    print *, ""
    print *, "Gallilean boosted: "
    print *, "t_prime_gal = t: ", t
    print *, "r_vec_prime_gal: ", [r_vec(1) + v_rel_vec(1)*t, r_vec(2) + v_rel_vec(2)*t, r_vec(3) + v_rel_vec(3)*t]
end subroutine test_lorentz_boost

program test
    use constants
    use user_variables

    implicit none
    integer :: num_tests_run, num_tests_failed

    num_tests_run = 0
    num_tests_failed = 0
    call test_spherical_cartesian_coord_change(num_tests_run, num_tests_failed)
    call test_lorentz_boost(num_tests_run, num_tests_failed)

    print *, " "
    print *, "TESTS FINISHED"
    print *, "-------------------------------------- "
    print *, "Summary:"
    print *, "Num tests run:      ", num_tests_run
    print *, "Num tests succeded: ", num_tests_run - num_tests_failed
    print *, "Num tests failed:   ", num_tests_failed
    print *, "-------------------------------------- "
end program test
