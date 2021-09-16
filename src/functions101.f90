!=============================================================================!
!=============================================================================!
!          file ~/Programme/SNR/Nuclei101/functions104.f90                    !
!=============================================================================!

!=============================================================================!
!=============================================================================!
double precision function R_p(r, t)
   use SNR_data
   implicit none
   double precision r, t, r_sh
   double precision t_shock

   r_sh = t_shock(t)                   ! either z or r coordinate
   if (r > r_sh) then
      R_p = 1.d0                       ! R=1
   else
      R_p = 4.d0                       ! R=4
   end if

end function R_p
!=============================================================================!
!=============================================================================!
!            diffusion coefficient : Bohm diffusion                           !
!=============================================================================!
!double precision function D_coef(En,t) ! yr
!  implicit none
!  double precision En,t     ! eV, yr
!  double precision R_L      ! yr
!  D_coef = R_L(En,t)/3.d0   ! yr
!end function D_coef
!=============================================================================!
!=============================================================================!
!                         Larmor Radius                                       !
!=============================================================================!
double precision function R_L(En, t) ! yr
   use SNR_data
   implicit none
   double precision En, t, B     ! eV, yr

   B = B0_turb                 ! constant in time

!  if (En>1.e11) B = B * En/1.e11
!  if (En>1.e11) B = B * sqrt(En/1.e11)

!  if (t<t_EDST) then
!     B = 100.d0*B
!  else
!     B = B/20.d0
!  end if

   if (B .eq. 0.d0) then
      R_L = 1.d100             ! yr
   else
      R_L = 3.523d-21*En/B     ! yr,             R_l/pc=1.08d-3*En/1.d18/(B/G)
   end if

end function R_L
!=============================================================================!
!=============================================================================!
double precision function dNdEdt(t)             ! injection rate
   use internal
   use SNR_data
   implicit none
   double precision t, v_snr, a, r_snr
   double precision v_shock, R_shock

   select case (inj_model)
   case (0)
      dNdEdt = 1.d0
      return
   case (1, 4, 5)
      a = 1.d0
   case (2)
      a = 3.d0
   end select

   v_snr = v_shock(t)
   r_snr = R_shock(t)
   dNdEdt = r_snr**2*v_snr**a

!  if (t<t_EDST) then
!     dNdEdt = t**2*v_snr**a
!  else
!     dNdEdt = t_EDST**(6.d0/5.d0)*t**(4.d0/5.d0)*v_snr**a
!  end if

end function dNdEdt
!=============================================================================!
!=============================================================================!
double precision function v_shock(t)   ! dimensionless
   use SNR_data
   implicit none
   double precision t, x, v
   double precision t_star

   select case (inj_model)
   case (0)
      v_shock = 3.d-2
   case (1, 2)
      x = t/t_EDST
      t_star = t/t_ch
      if (x < 1.d0) then
         v = 2.01d0*(1.d0 + 1.720*t_star**1.5d0)**(-5.d0/3.d0)
      else
         v = 2.d0/5.d0*1.42d0*(1.42d0*t_star - 0.254d0)**(-3.d0/5.d0)
      end if
      v = v*v_ch          ! pc/yr
      v_shock = v*3.262d0 ! yr
   case (4)
      v = 1.d-2*(t)**(-alpha_sh)
      v_shock = v*40.d0*6.d0
   case (5)
      v = 1.d-2*(1.56d0*t/t_EDST - 0.56d0)**(-alpha_sh)
      v_shock = v*4.d0
   case default
      call error('wrong case in v_shock', 0)
   end select

end function v_shock
!=============================================================================!
!=============================================================================!
double precision function t_shock(t)
   use internal
   use SNR_data
   implicit none
   double precision t, t_star, r, x
   double precision v_shock

   select case (inj_model)
   case (0)
      t_shock = v_shock(t)*t              ! planar shock moving to the right
   case (1, 2)
      x = t/t_EDST
      t_star = t/t_ch
      if (x < 1.d0) then
         r = 2.01d0*t_star*(1.d0 + 1.72d0*t_star**1.5d0)**(-2.d0/3.d0)
      else
         r = (1.42d0*t_star - 0.254d0)**(2/5.d0)
      end if
      r = r*R_ch          ! pc
      t_shock = r*3.262d0 ! yr
   case (4)
      t_shock = 1.d-2*(t)**(-alpha_sh + 1.d0)/((1.d0 - alpha_sh))
      t_shock = 40.d0*t_shock*6.d0
   case (5)
      t_shock = 1.d-2/1.56d0*t_EDST* &
                (1.56d0*t/t_EDST - 0.56d0)**(-alpha_sh + 1.d0)/((1.d0 - alpha_sh))
      t_shock = 4.d0*t_shock
   case default
      call error('wrong case in t_shock', 0)
   end select

end function t_shock
!=============================================================================!
!=============================================================================!
double precision function R_shock(t)     ! pc
   use SNR_data
   implicit none
   double precision t
   double precision x, t_star

   select case (inj_model)
   case (0)
      call error('inj_model 0 not in R_shock', 0)
   case (1, 2, 4)
!     call error('R_shock not needed ?!?',0)
!     R_shock = R_EDST*(t/t_EDST)**(2.d0/5.d0)
      x = t/t_EDST
      t_star = t/t_ch
      if (x < 1.d0) then
         R_shock = 2.01d0*t_star*(1.d0 + 1.72d0*t_star**1.5d0)**(-2.d0/3.d0)
      else
         R_shock = (1.42d0*t_star - 0.254d0)**(2.d0/5.d0)
      end if
      R_shock = R_shock*R_ch          ! pc
   case (5)
      R_shock = 1.d-2/1.56d0*t_EDST* &
                (1.56d0*t/t_EDST - 0.56d0)**(-alpha_sh + 1.d0)/((1.d0 - alpha_sh))
      R_shock = 4.d0*R_shock  ! pc
      R_shock = R_shock/3.262d0 ! yr
   case default
      call error('wrong case in R_shock', 0)
   end select

end function R_shock
!=============================================================================!
!=============================================================================!
double precision function tau_syn(m, E, t)     ! pc
   use particle_data
   use SNR_data
   implicit none
   double precision m, E, t, p_perp, B, chi, Psynch
   double precision, parameter :: B_cr = 4.14d13     !crit. B/Gauss, electrons

   p_perp = E**2 - m**2      ! we assume that p_perp = p
!  write(*,*) E,m
!  write(*,*) log10(E),log10(m)
!  write(*,*) (E/1.d9)**2,(m/1.d9)**2
   p_perp = sqrt(p_perp)

   B = B0_reg + B0_turb

   chi = p_perp/m*B/B_cr*(m_e/m)**2       ! dimensionless

!!  Psynch = dE/dt

! Classical value, from Jackson's
   Psynch = alpha_em*m**2*chi**2/1.5d0 ! eV^2
   Psynch = Psynch*4.8d22 ! eV/yr (1d7/197 eV/cm * 0.9461d18 cm/yr)

! new expression by Baier, V.N, Kathov, V.M, valid for all chi's

   Psynch = Psynch/(1.d0 + 4.8d0*(1.d0 + chi)*log(1.d0 + 1.7d0*chi) + 3.44d0*chi**2)**(2./3.d0)

   tau_syn = E/Psynch  ! eV/(eV/yr) = yr

end function tau_syn
!=============================================================================!
!=============================================================================!
double precision function cycle_energy_gain(theta_in, theta_out, v) result(e_gain)
   ! Computes the relative energy gain of a particle from a shock crossing cycle
   ! Computed in the rest frame of the SNR (v1 = 0)
   ! v = v2 - v1, function of v_shock
   ! Can use for comparison, but will probably be scrapped
   ! Assumed c = 1
   implicit none
   double precision, intent(in) :: theta_in, theta_out, v
   double precision :: beta

   e_gain = &
      (1 - v*cos(theta_in) + beta*cos(theta_out) - (beta**2)*cos(theta_in)*cos(theta_out)) &
      /(1 - v**2) - 1
end function cycle_energy_gain
!=============================================================================!
!=============================================================================!
!                   end file functions101.f90                                 !
!=============================================================================!

