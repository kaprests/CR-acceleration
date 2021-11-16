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
subroutine small_angle_dev(theta, phi, E, v_sh)
   ! Random small angle within a cone centered around z-axis
   use constants, only: pi, two_pi
   use particle_data, only: m_p
   use user_variables, only: theta_max ! Maximal scattering angle 
   implicit none
   double precision, intent(inout) :: phi, theta
   double precision, intent(in) :: E, v_sh
   double precision :: cos_theta_cone, theta_cone  ! Opening angle of loss cone
   double precision :: ran0

!   ! Compute loss cone opening, theta_cone
!   ! Set max scattering, theta_max, to 10% of loss cone opening
!   cos_theta_cone = v_sh/sqrt(1 - (E**2)/(m_p**2))
!   if (abs(cos_theta_cone) > 1) call error("cosine exceeds 1, small_angle", 0)
!   theta_cone = acos(cos_theta_cone)
!   theta_max = 0.1 * theta_cone

   ! Temporary max scattering values
   ! theta_max = pi => isotropic 
   ! theta_max = pi/10

   ! Random angle within the max scattering cone
   theta = theta_max * ran0() ! Theta within max
   phi = two_pi * ran0() ! Azimuthal angle phi isotropic
end subroutine small_angle_dev


subroutine euler_RyRz(theta, phi, R)
   implicit none
   double precision, intent(in) :: theta, phi
   double precision, intent(inout) :: R(3,3)
   double precision :: ct, cp, st, sp

   ct = cos(theta)
   cp = cos(phi)
   st = sin(theta)
   sp = sin(phi)

   ! R(column, row)
   R(1,1) = ct*cp
   R(1,2) = sp
   R(1,3) = -st*cp

   R(2,1) = -ct*sp
   R(2,2) = cp
   R(2,3) = st*sp

   R(3,1) = st
   R(3,2) = 0.d0
   R(3,3) = ct
end subroutine euler_RyRz


subroutine matrix_vec_mult(M, vec)
   implicit none
   double precision, intent(in) :: M
   double precision, intent(inout) :: vec
   integer :: i, j
   
   i = 0
   do i = 1, 3, 1
   do j = 1, 3, 1
      ! Look up or work it out
   end do
   end do
end subroutine matrix_vec_mult
