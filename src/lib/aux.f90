!=============================================================================!
!=============================================================================!
!                         error handling                                      !
!=============================================================================!
subroutine error(string, s)
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
!=============================================================================!
!=============================================================================!
!          random number generator from Numerical Recipes (Fortran90)         !
!=============================================================================!
function ran0()
   use internal, only: iseed
   implicit none
   integer, parameter :: K4B = selected_int_kind(9)
!  integer(K4B), intent(inout) :: iseed
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
!=============================================================================!
!=============================================================================!
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
!=============================================================================!
!=============================================================================!
subroutine scattering_angle_dev(theta, phi)
   ! Random small angle within a cone centered around z-axis
   use constants, only: pi, two_pi
   use user_variables, only: theta_max ! Maximal scattering angle
   implicit none
   double precision, intent(out) :: phi, theta
   double precision :: ran0, z

   ! Random angle within the max scattering cone
   z = cos(theta_max) + (1 - cos(theta_max)) * ran0()
   theta = acos(z) ! Theta within max
   phi = two_pi*ran0() ! Azimuthal angle phi isotropic
end subroutine scattering_angle_dev

subroutine scattering_angle(theta, phi, theta_max)
   ! Random small angle within a cone centered around z-axis
   use constants, only: pi, two_pi
   implicit none
   double precision, intent(inout) :: phi, theta
   double precision, intent(in) :: theta_max
   double precision :: ran0, z

   ! Random angle within the max scattering cone
   z = cos(theta_max) + (1 - cos(theta_max)) * ran0()
   theta = acos(z) ! Theta within max
   phi = two_pi*ran0() ! Azimuthal angle phi isotropic
end subroutine scattering_angle

subroutine max_scattering_angle(theta_max_computed, v_shock, E_particle)
   ! Computes the loss cone angle, and sets max_pitch scattering angle
   ! to some fraction of cone angle.
   use constants, only: pi
   use particle_data, only: m_p
   use user_variables, only: theta_max, theta_max_set

   implicit none
   double precision, intent(out) :: theta_max_computed
   double precision, intent(in) :: v_shock, E_particle
   double precision :: v_particle, v_p
   double precision :: cos_theta_cone, theta_cone

   if (theta_max_set) then
      ! Use user provided theta max
      theta_max_computed = theta_max
   else
      v_p = v_particle(E_particle, m_p)
      if (v_shock > v_p) then
         ! isotropic -- E.g. particles injected in front of UR shock (before overtaken 1st time)
         theta_max_computed = pi
      else 
         ! Compute loss cone opening, theta_cone
         ! Set max scattering, theta_max, to 100% of loss cone angle
         cos_theta_cone = v_shock / v_p
         if (abs(cos_theta_cone) > 1) then 
            print *, "v_shock: ", v_shock
            print *, "E: ", E_particle
            print *, "m: ", m_p
            print *, "cos_theta_cone: ", cos_theta_cone
            call error("cosine exceeds 1, max_scattering_angle", 0)
         endif
         theta_cone = acos(cos_theta_cone)
         theta_max_computed = 1.0*theta_cone
         !print *, "!!!!!!!!!!!!!!!!!"
         !print *, "theta max computed: ", theta_max_computed
         !print *, "!!!!!!!!!!!!!!!!!"
      endif
   endif
end subroutine max_scattering_angle

subroutine euler_RyRz(theta, phi, R)
   implicit none
   double precision, intent(in) :: theta, phi
   double precision, intent(inout) :: R(3, 3)
   double precision :: ct, cp, st, sp

   ct = cos(theta)
   cp = cos(phi)
   st = sin(theta)
   sp = sin(phi)

   ! R(column, row)
   R(1, 1) = ct*cp
   R(1, 2) = sp
   R(1, 3) = -st*cp

   R(2, 1) = -ct*sp
   R(2, 2) = cp
   R(2, 3) = st*sp

   R(3, 1) = st
   R(3, 2) = 0.d0
   R(3, 3) = ct
end subroutine euler_RyRz

subroutine radially_outward(phi_rad, theta_rad, x1, x2, x3)
   ! Finds the angles corresponding to radially out at point x, i.e. shock normal
   use constants, only: pi, two_pi
   implicit none
   double precision, intent(in) :: x1, x2, x3
   double precision, intent(inout) :: phi_rad, theta_rad
   theta_rad = atan2(sqrt(x1**2 + x2**2), x3)
   phi_rad = atan(x2/x1)
   if (x1 == 0) then
      if (x2 > 0) then
         phi_rad = pi/4
      else if (x2 < 0) then
         phi_rad = 3*pi/4
      else
         phi_rad = pi
      end if
   end if
   if (x1 < 0.d0 .and. x2 > 0) phi_rad = phi_rad + pi
   if (x1 < 0.d0 .and. x2 < 0) phi_rad = phi_rad + pi
   if (x1 > 0.d0 .and. x2 < 0) phi_rad = phi_rad + two_pi
end subroutine radially_outward

double precision function get_v_2(v_shock) result(v)
   use user_variables, only: gamma_sh
   implicit none
   double precision, intent(in) :: v_shock
   if (v_shock <= 0.1 .and. 0 < v_shock) then
      ! Non relativistic regime
      v = 0.75d0*v_shock
   else if (0.9 <= v_shock .and. v_shock < 1) then
      ! (Ultra) relativistic regime
      v = 1.d0 - 1.d0/(gamma_sh**2)
   else if (0.1 < v_shock .and. v_shock < 0.9) then
      ! In between regimes
      v = 0.75d0*v_shock
      call error("v_shock in between regimes 0.1 < v_shock < 0.9", 1)
   else
      call error("v_shock outside of allowed range [0, 1]", 0)
   end if
end function get_v_2

!subroutine matrix_vec_mult(M, vec)
!   implicit none
!   double precision, intent(in) :: M
!   double precision, intent(inout) :: vec
!   integer :: i, j
!
!   i = 0
!   do i = 1, 3, 1
!   do j = 1, 3, 1
!      ! Look up or work it out
!   end do
!   end do
!end subroutine matrix_vec_mult
