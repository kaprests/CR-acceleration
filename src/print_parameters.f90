module constants
    implicit none
    save
    double precision, parameter :: &
        N_av = 6.022142d23, & ! Avogadro's number
        pi = 3.1415926536d0, &
        two_pi = 2.d0*pi, &
        degree_rad = pi/180.d0, &
        rad_degree = 180.d0/pi, &
        E = 1.1d10, & ! Particle energy
        v_sh = 0.8, & ! Shock velocity (fraction of light speed)
        m_p = 938.272d6                            ! Proton mass
end module constants

program print_parameters
    use constants
    use particle_data, only: e_elm
    use SNR_data, only: B0_turb
    use stepsize_interpolated_polynom_coefficients
    use user_variables
    implicit none
    integer, parameter :: n = 3
    integer i, j
    double precision k(3), p(3), R_Euler(n, n), phi, theta, phi0, theta0
    double precision cubic_spline_small_angle_step_correction

  !!! Testing stuff !!!
    double precision :: R_L, D_coef, v_particle
    double precision :: E_inj, x_val
    double precision, dimension(100) :: x

    write(*, *) "theta_max: ", theta_max
!  print *, "!!!!!!!!!!!!!!!!!!!!!!!!"
!  print *, "E_inj(eV): ", E_inj
!  print *, "R_L(yr): ", R_L(E_inj, 1.0)
!  print *, "stepsize small angle corr: ", cubic_spline_small_angle_step_correction(pi*0.999)
!  print *, "stepsize (iso): ", R_L(E_inj, 1.0)*cubic_spline_small_angle_step_correction(pi*0.999)
!  print *, "R_L/3: ", R_L(E_inj, 1.0)/3
!  print *, "R_L: ", sqrt(E_inj**2 - m_p**2)/(abs(e_elm)*B0_turb)
!  print *, "R_L - R_L: ", sqrt(E_inj**2 - m_p**2)/(abs(e_elm)*B0_turb) - R_L(E_inj, 1.0)
!  print *, "D_coeff(?): ", D_coef(E_inj, 1.0)
!  print *, "v_particle: ", v_particle(E_inj, m_p)
!  print *, "!!!!!!!!!!!!!!!!!!!!!!!!"
!  print *, "Testing cubic spline stepsize: "
!  print *, "bp: ", bp
!  print *, " "
!  print *, "coeffs(:, 1): ", coeffs(:, 1)
!  print *, " "
!  print *, "coeffs(:, 2): ", coeffs(:, 2)
!  print *, " "
!  print *, "coeffs(:, 3): ", coeffs(:, 3)
!  print *, " "
!  print *, "coeffs(:, 4): ", coeffs(:, 4)
!  print *, "cs(0.0005pi): ", cubic_spline(0.005*pi)
!  print *, "cs(0.5pi): ", cubic_spline(0.5*pi)
!  print *, "cs(1.0pi): ", cubic_spline(0.99*pi)
!
!  print *, "///////////////"
!  x_val = pi/100
!  do i = 1, 100, 1
!     x(i) = i * pi/100
!     print *, x(i), ", ",cubic_spline_stepsize(x(i))
!     end do
end program print_parameters
