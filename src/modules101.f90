!==============================================================================!
!==============================================================================!
module user_variables
  ! Maybe read user_variables from file?
  implicit none
  save

  integer, parameter ::                      &
    n_sets = 10,                            & ! number of MC set
    n_start = 1*10**2,                       & ! injected particles/set
    debug = 0,                               & ! 0, no debugging info 
    restart = 0,                             & ! 1 use old data
    iseed_shift = 0                           ! positive shift of random seed 

  character(4) ::  basename = '_mod'            ! name in output
  character(len=:), allocatable :: filename     ! name in output
  character(10) :: outdir = './Data'

end module user_variables
!==============================================================================!
!==============================================================================!
module SNR_data
  implicit none
  save

  ! move to user variables?
  integer, parameter ::  inj_model = 2      ! stat. (0), SNR (1/2, Voelk/Russ) 

  double precision, parameter ::      &
    M_sun =  1.78698d54,              &  ! erg
    M_ej = 4.d0*M_sun,                &  ! ejected mass (erg)
    E_snr = 5.d51,                    &  ! mechanical explosion energy (erg)
    n_ISM = 2.d0,                     &  ! density of the ISM (1/cm^3)
    xi_scat = 1.d3,                   &  ! artif. enhanc. of interactions
    B0_reg = 1d-6,                    &  ! regular B/Gauss
    B0_turb = B0_reg/1.d0,            &  ! turbulent B/Gauss
    t_max = 1.3d4                        ! lifetime in yr \gsim t_inj_fin

  double precision, parameter ::      &
    alpha_sh =  0.6d0                    ! inj_model 4

  double precision t_ch,R_ch,v_ch,t_EDST,R_EDST !,v_EDST

end module SNR_data
!==============================================================================!
!==============================================================================!
module constants
  implicit none
  save
  double precision, parameter ::             & 
    N_av = 6.022142d23,                      & ! Avogadro's number
    pi = 3.1415926536d0,                     &
    two_pi = 2.d0*pi
end module constants
!==============================================================================!
!==============================================================================!
module event_internal                           ! used instead of stack
  implicit none
  save
  integer, parameter ::  n_max = 2              ! max. number of secondaries
  type one_event
     integer :: pid,A,Z
     double precision :: x(3),t,E,w
  end type one_event 
  type (one_event), target :: event(n_max)
end module event_internal
!=============================================================================!
!=============================================================================!
module stack                                      ! not used here
  implicit none
  save
  integer, parameter ::  n_maxs = 20000           ! max. number of secondaries

  type :: fix
     double precision :: x(3),r,z,t,w
  end type fix 

  type :: var
     integer :: pid,A,Z
     double precision :: E
  end type var
  
  type :: one_particle
     type(fix) :: fix
     type(var) :: var
  end type one_particle
  type (one_particle), target :: events(n_maxs)

end module stack
!==============================================================================!
!==============================================================================!
module particle
  implicit none
  save
  integer, parameter ::               &
    pid_max = 15,                     &  ! maximal pid 
    n_stable = 11                         ! number of stable particles
  integer stable_pid(n_stable)
   data stable_pid /-8,-7,-5,-4,-1,0,1,4,5,7,8/      !...gamma,e,nue,numu,p
end module particle
!==============================================================================!
!==============================================================================!
module particle_data
!  use particle
  implicit none
  save
  double precision, parameter ::             &   ! particle masses in eV
    m_p = 938.272d6,                  &   ! proton mass 
    m_pi = 139.57d6,                  &   ! charged pion mass
    m_n = 939.56536d6,                &   ! neutron mass 
    m_pi0 = 134.9766d6,               &   ! neutral pion mass
    m_mu = 105.658d6,                 &   ! muon mass 
    m_kaon_charged = 493.677d6,       &   ! charged kaon mass 
    m_kaon_0long = 497.648d6,         &   ! kaon0 mass
    m_eta = 547.85d6,                 &     ! eta mass
    m_Lambda = 1115d6,                &     ! Lambda mass
    m_e = 0.5109989d6,                &     ! electron mass
    Q = m_n - m_p
  double precision, parameter ::             &
    alpha_em = 7.2973526d-3,          &   ! alpha_em
    e_elm = 0.302822121d0                 ! e_em = sqrt(alpha*4*pi)
end module particle_data
!==============================================================================!
!==============================================================================!
module internal
  use particle
  implicit none
  save

  integer iseed,n_in,n_out

  double precision, parameter ::             & ! all energies in eV
    E_inj = 1.1d10,                          & ! initial energy   
    E_min = 1d10,                            & ! minimal energy for bining
    E_min_em = 1.d8                            ! minimal em energy to be stored

  double precision, parameter :: t_inj_init = 1.d0
  double precision t_inj_fin,d_f,d_em,d_time_out,d_time_inj,alpha_f, &
    alpha_f1,K_inj

  double precision, parameter ::  dn = 0.1d0

  integer n_tot

end module internal
!==============================================================================!
!==============================================================================!
module result
  use particle
  implicit none
  save

  integer, parameter ::                      &
    n_tbin_in = 200,                         & ! time bins for injected spectrum
    n_tbin_out = 20,                         & ! time bins for output spectrum
    n_enbin = 100                              ! ten decades starting in E_min

  double precision :: N_i(2,n_tbin_in),E2N_i(2,n_tbin_in),  &
    En_f(-pid_max:pid_max,n_enbin),            &
    En_f_tot(-pid_max:pid_max,n_enbin) 

end module result
!==============================================================================!
!==============================================================================!
module test_var
  implicit none
  save
  integer n_injected,sec,flag,accel
end module test_var
!==============================================================================!
!==============================================================================!
