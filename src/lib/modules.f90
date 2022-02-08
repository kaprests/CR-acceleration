!===================!
! File: modules.f90 !
!===================!


module constants
   implicit none
   save
   double precision, parameter :: &
      N_av = 6.022142d23, & ! Avogadro's number
      pi = 3.1415926536d0, &
      two_pi = 2.d0*pi
end module constants


module user_variables
   ! Maybe read user_variables from file? Add option to override with cmd-line flag
   use constants, only: pi
   implicit none
   save
   integer :: &
      n_sets = 10, & ! number of MC set
      n_start = 1*10**2, & ! injected particles/set
      debug = 0, & ! 0, no debugging info
      restart = 0, & ! 1 use old data
      iseed_shift = 0, &                        ! positive shift of random seed
      inj_model = 0                             ! constant shock_velocity
   double precision :: shock_velocity = 3.d-2   ! velocity of constant shock
   double precision :: gamma_sh                 ! gamma of constant velocity shock
   logical :: gamma_set = .false.               ! True if gamma was set instead of vshock (in init)

   character(4) ::  basename = '_pas'           ! name in output
   character(len=:), allocatable :: filename    ! name in output
   character(10) :: outdir = './Data'
   character(10) :: outdir_raw = './RawData'
   integer :: num_steps_log = 200               ! Must reduce for large particle count
   integer :: num_cross_log = 100               ! Must reduce for large particle count
   logical :: isotropic = .false.               ! Use isotropic rw if true
   double precision :: E_inj_exp ! initial energy exponent (base 10) (E_inj)
   character(10) :: E_inj_exp_str

   ! Strings
   character(10) :: n_start_str, n_sets_str, v_shock_str, gamma_str, n_proc_str

   ! For shockless random walks only
   integer :: num_steps_tot = 10000
   character(10) :: num_steps_tot_str
   double precision :: theta_max_pi_frac = 1.0  ! fraction of pi
   double precision :: theta_max
   character(10) :: theta_max_str
   logical :: theta_max_set = .false.
   double precision :: t_max = -1 ! -1 temp value, means not set (will default to t_max_str)
   character(10) :: t_max_str
   double precision :: stepsize_exp = 2.1
   character(10) :: stepsize_exp_str
   character(10) :: stepsize_str
end module user_variables


module SNR_data
   use user_variables, only: inj_model
   implicit none
   save
   !integer, parameter ::  inj_model = 2      ! stat. (0), SNR (1/2, Voelk/Russ)

   double precision, parameter :: &
      M_sun = 1.78698d54, &  ! erg
      M_ej = 4.d0*M_sun, &  ! ejected mass (erg)
      E_snr = 5.d51, &  ! mechanical explosion energy (erg)
      n_ISM = 2.d0, &  ! density of the ISM (1/cm^3)
      xi_scat = 1.d3, &  ! artif. enhanc. of interactions
      B0_reg = 1d-6, &  ! regular B/Gauss
      B0_turb = B0_reg/1.d0, & ! turbulent B/Gauss
      t_max_snr = 1.3d4 ! lifetime in yr \gsim t_inj_fin

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
!  use particle
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
   integer iseed, n_in, n_out
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
   use user_variables, only: n_sets, n_start

   implicit none
   save

   integer, parameter :: &
      n_tbin_in = 200, & ! time bins for injected spectrum
      n_tbin_out = 20, & ! time bins for output spectrum
      n_enbin = 100      ! ten decades starting in E_min

   double precision :: N_i(2, n_tbin_in), &
      E2N_i(2, n_tbin_in), &
      En_f(-pid_max:pid_max, n_enbin), &
      En_f_tot(-pid_max:pid_max, n_enbin), &
      ! Protons only
      NE_esc(n_enbin) = 0, & ! # protons escaped at each energy(bin)
      NE_esc_tot(n_enbin) = 0 ! # protons escaped at each energy(bin) (total for a proc)
   double precision :: rel_energy_gain_total_sum
   double precision, allocatable :: exit_energies(:) ! # unbinned exit energies
   double precision, allocatable :: num_crossings_total(:) ! # number of shock crossings
   double precision, allocatable :: trajectories(:, :, :)
   double precision, allocatable :: crossing_flight_angles(:, :)
   double precision, allocatable :: phase_space_dist(:, :, :)

   ! For shockless random walks only
   double precision, allocatable :: final_positions(:, :)
   double precision, allocatable :: sample_positions(:, :, :) ! Allocated in random walk/pitch angle routine
end module result


module test_var
   implicit none
   save
   integer n_injected, sec, flag, accel
end module test_var


module stepsize_interpolated_polynom_coefficients
   ! Consider reading from a file
   implicit none
   ! Cubic splines:
   integer, parameter :: N_bp = 18, N_coeffs = 4
   double precision, dimension(N_bp) :: bp = &
      [&
         0.01570796, &
         0.15707963, &
         0.31415927, &
         0.34557519, &
         0.37699112, &
         0.40840704, &
         0.43982297, &
         0.4712389, &
         0.62831853, &
         0.78539816, &
         0.9424778, &
         1.25663706, &
         1.57079633, &
         1.88495559, &
         2.19911486, &
         2.51327412, &
         2.82743339, &
         3.14159265 &
      ]

   double precision, dimension(N_bp-1, N_coeffs) :: coeffs = &
      reshape(&
         [&
            7.67837462e-07, &
            7.67837462e-07, &
            -1.79641382e-05, &
            3.23022279e-05, &
            -2.75310810e-05, &
            3.67902326e-05, &
            -3.26158472e-05, &
            1.73269275e-06, &
            -4.83699610e-08, &
            3.47486191e-06, &
            -2.47089611e-06, &
            2.99821719e-06, &
            -4.27385164e-06, &
            -3.36074962e-06, &
            -6.55789674e-06, &
            -8.27279487e-06, &
            -8.27279487e-06, &
            !-------------
            4.06638197e-06, &
            4.39203336e-06, &
            4.75386824e-06, &
            3.06078810e-06, &
            6.10520136e-06, &
            3.51045810e-06, &
            6.97785584e-06, &
            3.90388466e-06, &
            4.72039688e-06, &
            4.69760307e-06, &
            6.33509317e-06, &
            4.00632845e-06, &
            6.83208158e-06, &
            2.80407131e-06, &
            -3.63360592e-07, &
            -6.54403266e-06, &
            -1.43409581e-05, &
            !-------------
            1.66155357e-07, &
            1.36193565e-06, &
            2.79857051e-06, &
            3.04407518e-06, &
            3.33203324e-06, &
            3.63411809e-06, &
            3.96361819e-06, &
            4.30547815e-06, &
            5.66017712e-06, &
            7.13955310e-06, &
            8.87256497e-06, &
            1.21214184e-05, &
            1.55264053e-05, &
            1.85536920e-05, &
            1.93204639e-05, &
            1.71504423e-05, &
            1.05892290e-05, &
            !-------------
            1.08493499e-09, &
            1.08014581e-07, &
            4.33291985e-07, &
            5.25346550e-07, &
            6.25001441e-07, &
            7.34852308e-07, &
            8.53626906e-07, &
            9.84023216e-07, &
            1.76336618e-06, &
            2.76874837e-06, &
            4.01960328e-06, &
            7.35563712e-06, &
            1.16520653e-05, &
            1.70716126e-05, &
            2.30729733e-05, &
            2.89034778e-05, &
            3.33890695e-05 &
         ], shape(coeffs))
end module stepsize_interpolated_polynom_coefficients
