!=============================================================================!
!=============================================================================!
module constants
  implicit none
  save
  double precision, parameter ::             & 
    N_av = 6.022142d23,                      & ! Avogadro's number
    pi = 3.1415926536d0,                     &
    two_pi = 2.d0*pi,                        &
    degree_rad = pi/180.d0,                  &     
    rad_degree = 180.d0/pi,                  &
    theta_cone = 5.0                          ! Cone opening, degrees, temporary value

end module constants
!=============================================================================!
!=============================================================================!
program test
  use constants
  implicit none
  integer, parameter :: n=3
  integer i,j
  double precision k(3),p(3),R_Euler(n,n),phi,theta,phi0,theta0

  
! initial momentum vector
  theta0=80.d0*degree_rad
  phi0=00.d0*degree_rad 
  
  k(1) = cos(phi0)*sin(theta0)
  k(2) = sin(phi0)*sin(theta0)
  k(3) = cos(theta0)


! find scattering angles (theta,phi) in system where p=(1,0,0):
  call small_angle(theta,phi)

  write(*,*) 'scattering angles',theta*rad_degree,phi*rad_degree
  
  p(1) = cos(phi)*sin(theta)
  p(2) = sin(phi)*sin(theta)
  p(3) = cos(theta)

  
! and rotate p back to the old coordinates:
  call Euler_backward(-phi0,-theta0,R_Euler)
  k = 0.d0
  do i=1,3 ! for column i in columns
     do j=1,3 ! for row j in rows
        k(i) = k(i) + R_Euler(i,j)*p(j)
     end do
  end do

! Equivalent back rotation
!  k(1) = k(1)*cb*ca - k(2)*cb*sa + k(3)*sb
!  k(2) = k(2)*ca + k(1)*sa
!  k(3) = -k(1)*ca*sb + k(2)*sb*sa + k(3)*cb

!.. calculate new theta, phi from cartesian coordinates:
  theta = atan2( sqrt(k(1)**2+k(2)**2), k(3) )
  phi   = atan( k(2)/k(1) ) 
  if (k(1)<0.d0.and.k(2)>0) phi = phi+pi
  if (k(1)<0.d0.and.k(2)<0) phi = phi+pi
  if (k(1)>0.d0.and.k(2)<0) phi = phi+two_pi

  write(*,*) 'new angles',theta*rad_degree,phi*rad_degree

end program test
!=============================================================================!
!=============================================================================!
subroutine Euler_backward(alpha,beta,R)
  implicit none
  integer i,j,k
  double precision alpha,beta,R(3,3),ca,sa,cb,sb
      
  ca = cos(alpha)
  sa = sin(alpha)
  cb = cos(beta)
  sb = sin(beta)

  ! R(column, row)
  R(1,1) =  cb*ca
  R(1,2) =  sa
  R(1,3) = -sb*ca 

  R(2,1) = -cb*sa
  R(2,2) =  ca
  R(2,3) =  sb*sa

  R(3,1) =  sb
  R(3,2) =  0.d0
  R(3,3) =  cb

end subroutine Euler_backward
!=============================================================================!
!=============================================================================!
subroutine small_angle(theta,phi)
  ! Generates a random small angle in rotated frame where p=(1,0,0)
  ! I.e. an angle within some solid angle centered around the z-axis
  use constants, only: theta_cone, two_pi, degree_rad
  implicit none
  double precision phi,theta,r
  double precision ran0

  theta = theta_cone * degree_rad * ran0()
  phi = two_pi * ran0()
end subroutine small_angle
!=============================================================================!
!=============================================================================!
