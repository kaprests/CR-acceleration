!============================================================================!
!============================================================================!
subroutine init(myid)
  use user_variables, only : restart
  implicit none
  integer myid

  call init_general(myid)  
  call init_inject_spec

!  if (myid==0.and.restart>0) call read_old_results

end subroutine init
!============================================================================!
!============================================================================!
subroutine init_general(myid)
  use internal
  use user_variables
  use SNR_data
  implicit none
  integer myid
  double precision v_EDST

! initialisation for random number (NumRec):
  iseed = 15321 + 2*(1+iseed_shift)*(myid+1)            

  n_in = 0                         ! # of particles on stack
  n_out = 0                        ! # of escaped particles

  d_f = log10(E_min)-0.1d0         ! init internal variables
  d_em = log10(E_min_em)-0.1d0

  if (myid==0) then
!     call test_int
     open(unit=99,file='Data/error'//filename)
     write(*,*) 'iseed = ',iseed
  end if
! transition time between ED and ST phases (yr) (McKee)  ! init SNR_data variables
  t_ch = 423.d0/sqrt(E_snr/1.d51)*(M_ej/M_sun)**(5.d0/6.d0)/n_ISM**(1.d0/3.d0) ! yr
  R_ch = 3.07d0*(M_ej/M_sun)**(1.d0/3.d0)/n_ISM**(1.d0/3.d0)! pc
  v_ch = R_ch/t_ch


  t_EDST = 0.495d0*t_ch  ! yr
  R_EDST = 0.727d0*R_ch ! pc
  v_EDST = sqrt(2.d0*E_snr/M_ej) 

! start and end=transition ST-PDS, Truelove/McKee for zeta_m = 1.
  t_inj_fin = 1.33d4*(E_snr/1d51)**(3.d0/14.d0) /n_ISM**(4.d0/7.d0)          
!  t_inj_fin = t_max

  if (myid==0 .and. inj_model>0) then
     write(*,*)
     write(*,*)'t_EDST = ',real(t_EDST),'/yr, v_EDST = ',real(v_EDST)
     write(*,*)'v_ch = ',real(v_ch)
     write(*,*)'t_inj_fin = ',real(t_inj_fin),'/yr'
     write(*,*)'compared to'
     write(*,*)'t_max     =',real(t_max),'/yr'
  end if

end subroutine init_general
!=============================================================================!
!=============================================================================!
! acceptance-rejection method, comparing dN/dEdt with f(t)=K*t**alpha_f       !
! with  f(t_EDST) = dN/dEdt(t_EDST) => K =  dN/dEdt(t_EDST)/t_EDST**alpha_f   !
!=============================================================================!
subroutine init_inject_spec
  use SNR_data
  use internal; use result
  implicit none
  double precision :: dNdEdt

  t_inj_fin = t_max

  d_time_inj = log10(t_inj_fin/t_inj_init)/dble(n_tbin_in ) ! linear or log??
  d_time_out = log10(t_max/t_inj_init)/dble(n_tbin_out)     ! linear or log??

  select case (inj_model)
  case (0)                     ! stationary 
      alpha_f = 1.d0           ! does not matter
  case (1,4)                   ! thermal / Voelk injection 
     alpha_f = 0.9d0           ! close to 1
  case (2)                     ! pressure / Russian
     alpha_f = 0.d0
  end select
  alpha_f1 = alpha_f + 1.d0

  K_inj = dNdEdt(t_EDST)
  K_inj = K_inj/t_EDST**alpha_f
  K_inj = K_inj*3d0              ! so that f(t) > dN/dEdt(t)

end subroutine init_inject_spec
!=============================================================================!
!=============================================================================!
subroutine inject !(i)
  use internal
  use event_internal
  use SNR_data
  implicit none
  integer, parameter :: n=1
  double precision :: x(3),t,w,phi,theta
  double precision :: dNdEdt0,f0,dNdEdt
  double precision r,ran0,t_shock


  select case(inj_model)
  case(0)
     t = 1.d2            ! t_inj_init
     x(1) = 0.d0
     x(2) = 0.d0
     x(3) =  t_shock(t)
  case(1,2,4)    
     do 
        r = ran0() 
         t = t_inj_init**alpha_f1 + &
            r*(t_inj_fin**alpha_f1 -  t_inj_init**alpha_f1) ! yr
        t = t**(1.d0/alpha_f1)
        dNdEdt0 = dNdEdt(t)
        f0 = K_inj*t**alpha_f
        r = ran0()
        if (dNdEdt0 .gt. f0) call error('dNdEdt0 > f0',0)
        if (dNdEdt0 .ge. r*f0) exit
     end do
     call isotropic(phi,theta)
     r = t_shock(t)                     ! delta function at shock
     x(1) = r*cos(phi)*sin(theta)
     x(2) = r*sin(phi)*sin(theta)
     x(3) = r*cos(theta)
!!!!!
!     t = t_inj_init
!     t = 2.d3
     x(1) = 0.d0
     x(2) = 0.d0
     x(3) = t_shock(t)               ! delta function at shock, planar approx
!!!!!
  case default
     write(*,*) 'injection model',inj_model
     call error('wrong shock model',0)
  end select

  w = 1.d0

  event(n)%pid = 7
  event(n)%A = 1
  event(n)%Z = 1 
  event(n)%x = x
  event(n)%E = E_inj
  event(n)%t = t
  event(n)%w = w


end subroutine inject
!=============================================================================!
!=============================================================================!
