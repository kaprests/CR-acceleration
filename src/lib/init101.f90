!============================================================================!
!============================================================================!
subroutine init(myid, n_proc)
   use user_variables
   implicit none
   integer, intent(in) :: myid, n_proc

   call init_general(myid, n_proc)
   call init_inject_spec

! if (myid==0.and.restart>0) call read_old_results

end subroutine init
!============================================================================!
!============================================================================!
subroutine parse_cmd_arguments
   ! Very naive, non flexible argument parser
   use user_variables
   use constants, only: pi
   implicit none

   integer :: i, n_args
   integer :: j, n_flags
   character(20) :: flag
   character(20) :: arg
   character(20), dimension(14), parameter :: flags = &
                                              [ &
                                              ! Integer:
                                              '--nsets     ', &  ! j=1
                                              '--nstart    ', &  ! j=2
                                              '--debug     ', &  ! j=3
                                              '--restart   ', &  ! j=4
                                              '--iseedshift', &  ! j=5
                                              '--injmod    ', &  ! j=6
                                              ! Float (double precision):
                                              '--vshock    ', &  ! j=7
                                              '--gamma     ', &  ! j=8
                                              ! Character
                                              '--fname     ', &  ! j=9
                                              ! Integer (shockless random walks)
                                              '--nsteps    ', &  ! j=10
                                              ! Float (shockless random walks)
                                              '--max-pi-fr ', &  ! j=11
                                              '--t-max     ', &  ! j=12
                                              '--iso       ', &  ! j=13
                                              '--stepexp   ' &  ! j=14
                                              ]

   n_args = command_argument_count()
   n_flags = size(flags)
   if (modulo(n_args, 2) .ne. 0) then
      ! Sort of crude error handling
      call error('Argument error, odd number of provided arguments and flags', 0)
   end if
   print *, 'Provided parameters/settings: '
   do i = 1, n_args, 2
      call get_command_argument(i, flag)
      call get_command_argument(i + 1, arg)
      if (flag(1:1) == '#') then
         ! Comment line, ignore
         cycle
      end if
      print *, trim(flag)//' '//trim(arg)
      do j = 1, n_flags, 1
         if (flag == trim(flags(j))) then
            ! integer args
            select case (j)
            case (1)
               read (arg, *) n_sets
            case (2)
               read (arg, *) n_start
            case (3)
               read (arg, *) debug
            case (4)
               read (arg, *) restart
            case (5)
               read (arg, *) iseed_shift
            case (6)
               read (arg, *) inj_model
            case (7)
               read (arg, *) shock_velocity
               gamma_set = .false.
            case (8)
               read (arg, *) gamma_sh
               shock_velocity = sqrt(1 - 1/(gamma_sh**2))
               !shock_velocity = 1 - 1/(2*gamma_sh**2)
               gamma_set = .true.
            case (9)
               basename = arg
            case (10)
               read (arg, *) num_steps_tot
            case (11)
               read (arg, *) theta_max_pi_frac
               theta_max = pi*theta_max_pi_frac
               theta_max_set = .true.
            case (12)
               read (arg, *) t_max
            case (13)
               if (arg == "true") then
                  isotropic = .true.
               else if (arg == "false") then
                  isotropic = .false.
               else
                  write (*, *) "Unrecognized argument for the --iso flag. Using default (false)."
               end if
            case (14)
               read (arg, *) stepsize_exp
            end select
            exit
         elseif (j == n_flags) then
            call error('Argument error, invalid flag: '//flag, 0)
         end if
      end do
   end do
end subroutine parse_cmd_arguments
!============================================================================!
!============================================================================!
subroutine init_general(myid, n_proc)
   use internal
   use user_variables
   use SNR_data
   use result
   implicit none
   integer, intent(in) :: myid, n_proc
   double precision v_EDST, v_shock
   !character(10) :: n_start_str, n_sets_str, , v_shock_str, gamma_str

   ! Parse command line arguments, and apply given settings/config
   ! Default values in are code overridden by values in file (future)
   ! Values set in file are overridden by command line arguments
   ! TODO: read parameters from file
   call parse_cmd_arguments ! Command line arguments

   ! Adds configuration metadata to filename
   write (n_sets_str, '(I10)') n_sets
   write (n_start_str, '(I10)') n_start
   write (n_proc_str, '(I10)') n_proc
   n_sets_str = adjustl(n_sets_str)
   n_start_str = adjustl(n_start_str)
   n_proc_str = adjustl(n_proc_str)
   if (isotropic) then
      basename = '_iso'
   end if
   filename = trim(basename)
   if (inj_model == 0) then
      ! constant shock velocity
      if (gamma_set) then
         ! gamma instead of vshock
         write (gamma_str, '(f10.3)') gamma_sh
         filename = filename//'_gamma'//trim(adjustl(gamma_str))
      else
         print *, v_shock(0)
         write (v_shock_str, '(f10.3)') v_shock(0)
         filename = filename//'_vshock'//trim(adjustl(v_shock_str))
      end if
   else if (inj_model == 1) then
      filename = filename//'_injmod1'
   else if (inj_model == 2) then
      filename = filename//'_injmod2'
   end if
   if (theta_max_set) then
      write (theta_max_str, '(f10.3)') pi*theta_max_pi_frac
      theta_max_str = adjustl(theta_max_str)
      filename = filename//'_theta-max'//trim(theta_max_str)
   end if
   filename = filename//'_nsets'//trim(n_sets_str)//'_nstart'//trim(n_start_str)
   filename = filename//'_n_proc'//trim(n_proc_str)
   print *, "filename metadata: ", filename

   ! If no t_max given by user, set default exit time t_max_snr from SNR_data
   if (t_max == -1) then
      t_max = t_max_snr
   end if

   ! For shockless random walk
   write (num_steps_tot_str, '(I10)') num_steps_tot
   num_steps_tot_str = adjustl(num_steps_tot_str)
   write (t_max_str, '(f10.3)') t_max
   t_max_str = adjustl(t_max_str)
   write (stepsize_exp_str, '(f10.3)') stepsize_exp
   stepsize_exp_str = adjustl(stepsize_exp_str)
   print *, "=========================="
   print *, "For shockless random walk:"
   print *, "t_max: ", t_max
   print *, "theta max: ", pi*theta_max_pi_frac
   print *, "stepsize exp: ", stepsize_exp
   print *, "=========================="

   ! Allocate dynamic arrays
   ! Not finalized
   allocate (exit_energies(n_sets*n_start))
   allocate (num_crossings_total(n_sets*n_start))

   allocate (trajectories(4, n_start, num_steps_log))
   allocate (crossing_flight_angles(n_start, num_cross_log))

   ! For shockless random walks only
   allocate (final_positions(3, n_start))
   print *, "!!!!!!!!!"
   print *, "final_positions allocated"
   print *, "size final_positions(bytes): ", sizeof(final_positions)
   print *, "count final_positions(bytes): ", size(final_positions)
   print *, "n_start * 3: ", n_start * 3
   print *, "n_start : ", n_start
   print *, "!!!!!!!!!"

   ! initialisation for random number (NumRec):
   iseed = 15321 + 2*(1 + iseed_shift)*(myid + 1)

   n_in = 0                         ! # of particles on stack
   n_out = 0                        ! # of escaped particles

   d_f = log10(E_min) - 0.1d0         ! init internal variables
   d_em = log10(E_min_em) - 0.1d0

   if (myid == 0) then
      !   call test_int
      open (unit=99, file=trim(outdir)//'/error'//filename)
      write (*, *) 'iseed = ', iseed
   end if
   ! transition time between ED and ST phases (yr) (McKee)  ! init SNR_data variables
   t_ch = 423.d0/sqrt(E_snr/1.d51)*(M_ej/M_sun)**(5.d0/6.d0)/n_ISM**(1.d0/3.d0) ! yr
   R_ch = 3.07d0*(M_ej/M_sun)**(1.d0/3.d0)/n_ISM**(1.d0/3.d0)! pc
   v_ch = R_ch/t_ch

   t_EDST = 0.495d0*t_ch  ! yr
   R_EDST = 0.727d0*R_ch ! pc
   v_EDST = sqrt(2.d0*E_snr/M_ej)

   ! start and end=transition ST-PDS, Truelove/McKee for zeta_m = 1.
   t_inj_fin = 1.33d4*(E_snr/1d51)**(3.d0/14.d0)/n_ISM**(4.d0/7.d0)
   !  t_inj_fin = t_max

   if (myid == 0 .and. inj_model > 0) then
      write (*, *)
      write (*, *) 't_EDST = ', real(t_EDST), '/yr, v_EDST = ', real(v_EDST)
      write (*, *) 'v_ch = ', real(v_ch)
      write (*, *) 't_inj_fin = ', real(t_inj_fin), '/yr'
      write (*, *) 'compared to'
      write (*, *) 't_max     =', real(t_max), '/yr'
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
   use user_variables, only: t_max
   implicit none
   double precision :: dNdEdt

   t_inj_fin = t_max

   d_time_inj = log10(t_inj_fin/t_inj_init)/dble(n_tbin_in) ! linear or log??
   d_time_out = log10(t_max/t_inj_init)/dble(n_tbin_out)     ! linear or log??

   select case (inj_model)
   case (0)                     ! stationary
   case (1, 4)                   ! thermal / Voelk injection
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
   integer, parameter :: n = 1
   double precision :: x(3), t, w, phi, theta
   double precision :: dNdEdt0, f0, dNdEdt
   double precision r, ran0, t_shock

   select case (inj_model)
   case (0)
      t = 1.d2            ! t_inj_init
      x(1) = 0.d0
      x(2) = 0.d0
      x(3) = t_shock(t)*1.01
   case (1, 2, 4)
      do
         r = ran0()
         t = t_inj_init**alpha_f1 + &
             r*(t_inj_fin**alpha_f1 - t_inj_init**alpha_f1) ! yr
         t = t**(1.d0/alpha_f1)
         dNdEdt0 = dNdEdt(t)
         f0 = K_inj*t**alpha_f
         r = ran0()
         if (dNdEdt0 .gt. f0) call error('dNdEdt0 > f0', 0)
         if (dNdEdt0 .ge. r*f0) exit
      end do
      call isotropic(phi, theta)
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
      write (*, *) 'injection model', inj_model
      call error('wrong shock model', 0)
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
