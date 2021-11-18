module rel_acceleration
   implicit none
   private

   public small_angle_random_walk
contains


   subroutine small_angle_random_walk(set, n_injected)
      ! Small angle random walk
      use user_variables, only: debug, t_max
      use constants;
      use particle_data, only: m_p
      use event_internal
      use result
      use internal
      use test_var, only: sec, accel

      implicit none
      
      ! Variables and parameters
      integer, intent(in) :: set, n_injected
      integer :: k, n_step
      double precision :: r1, r_sh1, r2, r_sh2
      double precision :: phi_p, theta_p, phi_sh, theta_sh
      double precision :: gamma_sh, cos_theta
      integer :: num_crossings
      double precision, pointer :: E, x(:), t, w
      integer, pointer :: pid, A, Z
      
      ! Functions
      double precision :: ran0, R_L, t_shock, v_shock

      pid => event(n_in)%pid
      A => event(n_in)%A
      Z => event(n_in)%Z
      E => event(n_in)%E
      x => event(n_in)%x
      t => event(n_in)%t
      w => event(n_in)%w

      d1 = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
      if (sec == 0 .and. abs(d1/t_shock(t) - 1.01d0) .gt. 1.d-6) then
         call error('wrong initial condition, shock', 0)
      end if

      ! From old sim
      r = ran0()
      m = A*m_p
      f = 0.d0

      rel_energy_gain_sum = 0
      num_crossings = 0

      do
         ! Step size (also num steps?)

         ! Flight direction

         ! Check shock cross

         ! Exit
      end do
   end subroutine small_angle_random_walk


end module rel_acceleration
