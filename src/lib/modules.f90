! File: modules.f90

module constants
    implicit none
    save
    double precision, parameter :: &
        N_av = 6.022142d23, & ! Avogadro's number
        pi = 3.1415926536d0, &
        two_pi = 2.d0*pi
end module constants

module user_variables
    use constants, only: pi

    implicit none
    save

    integer :: &
        n_sets = 10, & ! number of MC set
        n_start = 100, & ! injected particles/set
        debug = 0, & ! 0, no debugging info
        restart = 0, & ! 1 use old data
        iseed_shift = 0, & ! positive shift of random seed,
        inj_model = 1, & ! 0: constant vel, 1, 2
        num_sample_pos_target = 500, & ! May have to reduce for large n_start
        num_traj_pos = 500 ! May also have to reduce for large n_start

    ! inj_model = 0: constant shock velocity
    double precision :: v_shock_const = 3.d-2 ! Default 3.d-2
    double precision :: gamma_shock
    logical :: gamma_shock_set = .false. ! if inj_model = 0, use v_shock_const by default

    ! Shockless simulation
    logical :: shockless = .false.
    double precision :: shockless_t_max = 120

    ! Initial step in z-direction
    logical :: init_z = .false.

    ! Injection energy xponent
    double precision :: E_inj_exp

    ! Max scattering angle
    double precision :: theta_max_pi_fraction = 1.0 ! fraction of pi theta/pi, default 1
    double precision :: theta_max = pi ! pi as default i.e. isotropic rw

    ! Turn on/off small angle stepsize correction
    logical :: no_stepsize_corr

    ! Filename
    character(len=6), parameter ::  basename = '_accel'            ! name in output
    character(len=6), parameter :: basename_rw = '_randw'
    character(len=6), parameter :: default_outdir = './Data'
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: outdir
    character(len=10) :: n_start_str, n_sets_str, v_shock_str, gamma_shock_str, n_proc_str
    character(len=10) :: E_inj_exp_str
    character(len=10) :: theta_max_str
    character(len=10) :: shockless_t_max_str
end module user_variables

module SNR_data
    use user_variables, only: inj_model
    implicit none
    save

    double precision, parameter :: &
        M_sun = 1.78698d54, &  ! erg
        M_ej = 4.d0*M_sun, &  ! ejected mass (erg)
        E_snr = 5.d51, &  ! mechanical explosion energy (erg)
        n_ISM = 2.d0, &  ! density of the ISM (1/cm^3)
        xi_scat = 1.d3, &  ! artif. enhanc. of interactions
        B0_reg = 1d-6, &  ! regular B/Gauss
        B0_turb = B0_reg/1.d0, &  ! turbulent B/Gauss
        t_max = 1.3d4                        ! lifetime in yr \gsim t_inj_fin
    double precision, parameter :: &
        alpha_sh = 0.6d0                    ! inj_model 4
    double precision t_ch, R_ch, v_ch, t_EDST, R_EDST !,v_EDST
end module SNR_data

module event_internal                           ! used instead of stack
    implicit none
    save
    integer, parameter ::  n_max = 2              ! max. number of secondaries
    type one_event
        integer :: pid, A, Z
        double precision :: x(3), t, E, w
    end type one_event
    type(one_event), target :: event(n_max)
end module event_internal

module stack                                      ! not used here
    implicit none
    save
    integer, parameter ::  n_maxs = 20000           ! max. number of secondaries

    type :: fix
        double precision :: x(3), r, z, t, w
    end type fix

    type :: var
        integer :: pid, A, Z
        double precision :: E
    end type var

    type :: one_particle
        type(fix) :: fix
        type(var) :: var
    end type one_particle
    type(one_particle), target :: events(n_maxs)
end module stack

module particle
    implicit none
    save
    integer, parameter :: &
        pid_max = 15, &  ! maximal pid
        n_stable = 11                         ! number of stable particles
    integer stable_pid(n_stable)
    data stable_pid/-8, -7, -5, -4, -1, 0, 1, 4, 5, 7, 8/      !...gamma,e,nue,numu,p
end module particle

module particle_data
    !use particle
    implicit none
    save
    double precision, parameter :: &   ! particle masses in eV
        m_p = 938.272d6, &   ! proton mass
        m_pi = 139.57d6, &   ! charged pion mass
        m_n = 939.56536d6, &   ! neutron mass
        m_pi0 = 134.9766d6, &   ! neutral pion mass
        m_mu = 105.658d6, &   ! muon mass
        m_kaon_charged = 493.677d6, &   ! charged kaon mass
        m_kaon_0long = 497.648d6, &   ! kaon0 mass
        m_eta = 547.85d6, &     ! eta mass
        m_Lambda = 1115d6, &     ! Lambda mass
        m_e = 0.5109989d6, &     ! electron mass
        Q = m_n - m_p
    double precision, parameter :: &
        alpha_em = 7.2973526d-3, &   ! alpha_em
        e_elm = 0.302822121d0                 ! e_em = sqrt(alpha*4*pi)
end module particle_data

module internal
    use particle

    implicit none
    save

    integer :: iseed, n_in, n_out
    double precision :: E_inj = 1.0d10 ! initial energy
    double precision, parameter :: & ! all energies in eV
        E_min = 1d10, & ! minimal energy for bining
        E_min_em = 1.d8                            ! minimal em energy to be stored

    double precision, parameter :: t_inj_init = 1.d0
    double precision t_inj_fin, d_f, d_em, d_time_out, d_time_inj, alpha_f, &
        alpha_f1, K_inj
    double precision, parameter ::  dn = 0.1d0
    integer n_tot
end module internal

module result
    use particle

    implicit none
    save

    integer, parameter :: &
        n_tbin_in = 200, & ! time bins for injected spectrum
        n_tbin_out = 20, & ! time bins for output spectrum
        n_enbin = 100                              ! ten decades starting in E_min
    double precision :: N_i(2, n_tbin_in), E2N_i(2, n_tbin_in), &
        En_f(-pid_max:pid_max, n_enbin), &
        En_f_tot(-pid_max:pid_max, n_enbin)

    ! Currently only for protons
    integer, parameter :: n_angle_bins = 100
    double precision, dimension(n_angle_bins) :: cross_angle_distribution

    ! For shockless simulation only:
    double precision, dimension(:, :), allocatable :: final_positions
    double precision, dimension(:, :, :), allocatable :: sample_positions, initial_trajectories
end module result

module test_var
    implicit none
    save
    integer n_injected, sec, flag, accel
end module test_var

module stepsize_powerlaw_params
    implicit none
    double precision, parameter :: a = 0.12428921948582322
    double precision, parameter :: b = 1.9975134238662615
end module stepsize_powerlaw_params

module stepsize_interpolated_polynom_coefficients
    ! Consider reading from a file
    implicit none
    ! Cubic splines:
    integer, parameter :: N_bp = 20, N_coeffs = 4
    double precision, dimension(N_bp) :: bp = &
        [ &
        0.0, &
        0.03141593, &
        0.06283185, &
        0.09424778, &
        0.12566371, &
        0.15707963, &
        0.18849556, &
        0.21991149, &
        0.25132741, &
        0.28274334, &
        0.31415927, &
        0.62831853, &
        0.9424778, &
        1.25663706, &
        1.57079633, &
        1.88495559, &
        2.19911486, &
        2.51327412, &
        2.82743339, &
        3.14159265 &
        ]
    double precision, dimension(N_bp - 1, N_coeffs) :: coeffs = &
        reshape( &
        [&
        2.70636020e-02, &
        2.70636020e-02, &
        -2.37073365e-02, &
        4.39546824e-02, &
        -1.87818231e-02, &
        -1.60417636e-01, &
        5.29602274e-01, &
        -1.10605100e+00, &
        1.51325172e+00, &
        -8.55627356e-01, &
        2.13457194e-02, &
        2.52729654e-03, &
        3.31140405e-02, &
        -5.01524665e-02, &
        1.83658893e-03, &
        -1.22393254e-01, &
        -3.07919793e-01, &
        -6.14881702e-02, &
        -6.14881702e-02, &
        !---------------
        1.20934365e-01, &
        1.23485050e-01, &
        1.26035734e-01, &
        1.23801370e-01, &
        1.27944002e-01, &
        1.26173856e-01, &
        1.11054850e-01, &
        1.60968689e-01, &
        5.67258375e-02, &
        1.99346452e-01, &
        1.18705474e-01, &
        1.38823340e-01, &
        1.41205261e-01, &
        1.72414509e-01, &
        1.25146923e-01, &
        1.26877867e-01, &
        1.15249434e-02, &
        -2.78682625e-01, &
        -3.36633860e-01, &
        !---------------
        1.14817364e-04, &
        7.79347976e-03, &
        1.56324064e-02, &
        2.34812705e-02, &
        3.13900846e-02, &
        3.93734326e-02, &
        4.68261922e-02, &
        5.53720637e-02, &
        6.22111390e-02, &
        7.02558872e-02, &
        8.02477831e-02, &
        1.61152846e-01, &
        2.49126426e-01, &
        3.47652982e-01, &
        4.41134663e-01, &
        5.20310586e-01, &
        5.63791112e-01, &
        4.79861051e-01, &
        2.86553676e-01, &
        !---------------
        0.00000000e+00, &
        1.23803670e-04, &
        4.91357058e-04, &
        1.10612080e-03, &
        1.96735659e-03, &
        3.07919850e-03, &
        4.43570601e-03, &
        6.03282196e-03, &
        7.89696186e-03, &
        9.95428889e-03, &
        1.23316599e-02, &
        4.99198565e-02, &
        1.14327193e-01, &
        2.07555711e-01, &
        3.32235706e-01, &
        4.83230700e-01, &
        6.55418476e-01, &
        8.24128697e-01, &
        9.45470100e-01 &
        ], shape(coeffs))
end module stepsize_interpolated_polynom_coefficients
