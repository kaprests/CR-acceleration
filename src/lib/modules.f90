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
    double precision :: gamma_shock_const
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
        double precision, dimension(3) :: x, p ! 3-position and 3-momentum
        double precision :: t, E, w 
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
    double precision :: &
        N_i(2, n_tbin_in), &
        E2N_i(2, n_tbin_in), &
        En_f(-pid_max:pid_max, n_enbin), &
        En_f_tot(-pid_max:pid_max, n_enbin)
    double precision, dimension(n_enbin) :: & ! Only for protons
        exit_energy_enumerate_dist, &
        exit_energy_enumerate_dist_tot

    ! Currently only for protons
    integer, parameter :: n_angle_bins = 100
    double precision, dimension(n_angle_bins) :: &
        cross_angle_distribution_first = 0.d0, &            ! The first shock cross
        cross_angle_distribution_first_tot = 0.d0, &
        cross_angle_distribution_updown = 0.d0, &           ! upstream -> downstream
        cross_angle_distribution_updown_tot = 0.d0, &
        cross_angle_distribution_downup = 0.d0, &           ! downstream -> upstream
        cross_angle_distribution_downup_tot = 0.d0

    ! To be implemented into simulation:
    integer, parameter :: n_angle_bins_aniso = 3000 ! fine resolution to capture anisotropy
    double precision, dimension(n_angle_bins_aniso) :: &
        cross_angle_distribution_aniso_updown = 0.d0, & ! us -> ds in small cone 
        cross_angle_distribution_aniso_updown_tot = 0.d0 ! to capture anisotropic dist.

    ! May be useful with planar shocks -- implement that first!
    double precision, dimension(n_angle_bins) :: &
        flight_angle_distribution_first = 0.d0, &
        flight_angle_distribution_first_tot = 0.d0, &
        flight_angle_distribution_updown = 0.d0, &
        flight_angle_distribution_updown_tot = 0.d0, &
        flight_angle_distribution_downup = 0.d0, &
        flight_angle_distribution_downup_tot = 0.d0

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
        -8.52953258e-02, &
        -8.52953258e-02, &
        3.17431882e-01, &
        -5.23209231e-01, &
        6.18452040e-01, &
        -5.83329140e-01, &
        5.01806319e-01, &
        -5.05165162e-01, &
        5.89587185e-01, &
        -3.62003051e-01, &
        1.58993899e-02, &
        1.93515405e-02, &
        -1.38188633e-02, &
        3.45772912e-02, &
        -9.97533617e-02, &
        -3.60455252e-02, &
        -3.75974150e-01, &
        -3.61723259e-02, &
        -3.61723259e-02, &
        !-------------
        1.32448346e-01, &
        1.24409451e-01, &
        1.16370556e-01, &
        1.46287806e-01, &
        9.69764975e-02, &
        1.55264229e-01, &
        1.00286753e-01, &
        1.47580884e-01, &
        9.99701892e-02, &
        1.55537472e-01, &
        1.21419489e-01, &
        1.36404311e-01, &
        1.54642708e-01, &
        1.41618736e-01, &
        1.74207065e-01, &
        8.01917366e-02, &
        4.62196295e-02, &
        -3.08127659e-01, &
        -3.42219273e-01, &
        !-------------
        -1.25407502e-04, &
        7.94401816e-03, &
        1.55083452e-02, &
        2.37600009e-02, &
        3.14023744e-02, &
        3.93267505e-02, &
        4.73551214e-02, &
        5.51421129e-02, &
        6.29191592e-02, &
        7.09461691e-02, &
        7.96470287e-02, &
        1.60644764e-01, &
        2.52079881e-01, &
        3.45153159e-01, &
        4.44372761e-01, &
        5.24294501e-01, &
        5.64007803e-01, &
        4.81726969e-01, &
        2.77414455e-01, &
        !-------------
        0.00000000e+00, &
        1.24136794e-04, &
        4.93848001e-04, &
        1.10575255e-03, &
        1.98035250e-03, &
        3.08177505e-03, &
        4.45241414e-03, &
        6.05465736e-03, &
        7.91699113e-03, &
        1.00106023e-02, &
        1.23817269e-02, &
        4.98801830e-02, &
        1.14410809e-01, &
        2.08438191e-01, &
        3.31920576e-01, &
        4.85624964e-01, &
        6.57133909e-01, &
        8.27226322e-01, &
        9.47032763e-01 &
        ], shape(coeffs))
end module stepsize_interpolated_polynom_coefficients
