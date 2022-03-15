! File: init.f90

subroutine init(myid, n_proc)
    !use user_variables, only: restart

    implicit none
    integer, intent(in) :: myid, n_proc

    call init_general(myid, n_proc)
    call init_inject_spec
    !if (myid==0.and.restart>0) call read_old_results
end subroutine init

subroutine parse_cmd_args
    ! naive and non flexible argument parser
    use user_variables
    use constants, only: pi
    use internal, only: E_inj

    implicit none
    integer :: i, j, n_args
    character(len=20) :: flag
    character(len=20) :: arg
    integer, parameter :: n_flags = 15
    character(len=20), dimension(n_flags), parameter :: flags = [ &
                                                        '--nsets         ', & ! j=1
                                                        '--nstart        ', & ! j=2
                                                        '--debug         ', & ! j=3
                                                        '--restart       ', & ! j=4
                                                        '--iseed_shift   ', & ! j=5
                                                        '--injmod        ', & ! j=6
                                                        '--vshock        ', & ! j=7
                                                        '--gamma         ', & ! j=8
                                                        '--fname         ', & ! j=9
                                                        '--max-pi-fr     ', & ! j=10
                                                        '--t-max         ', & ! j=11
                                                        '--E-inj-exp     ', & ! j=12
                                                        '--shockless     ', & ! j=13
                                                        '--init-z-ax     ', & ! j=14
                                                        '--no-step-corr  ' & ! j=15
                                                        ]

    n_args = command_argument_count()
    if (modulo(n_args, 2) .ne. 0) then
        call error('Argument error, odd number of provided arguments and flags', 0)
    end if

    write (*, *) "Provided parameters/settings:"
    do i = 1, n_args, 2
        call get_command_argument(i, flag)
        call get_command_argument(i + 1, arg)
        if (flag(1:1) == '#') then
            ! Comment line, skip
            cycle
        end if
        write (*, *) trim(flag)//' '//trim(arg)
        do j = 1, n_flags, 1
            if (flag == trim(flags(j))) then
                select case (j)
                case (1) ! n_sets
                    read (arg, *) n_sets
                case (2) ! n_start
                    read (arg, *) n_start
                case (3) ! debug
                    read (arg, *) debug
                case (4) ! restart
                    read (arg, *) restart
                case (5) ! iseed_shift
                    read (arg, *) iseed_shift
                case (6) ! inj_model
                    read (arg, *) inj_model
                case (7) ! v_shock_const
                    read (arg, *) v_shock_const
                    gamma_shock_set = .false.
                    inj_model = 0
                case (8) ! gamma_shock
                    read (arg, *) gamma_shock
                    gamma_shock_set = .true.
                    inj_model = 0
                case (9) ! fname
                    read (arg, *) filename
                case (10) ! max-pi-fr
                    read (arg, *) theta_max_pi_fraction
                    theta_max = pi*theta_max_pi_fraction
                case (11) ! t-max
                    read (arg, *) shockless_t_max
                case (12) ! E_inj_exp
                    read (arg, *) E_inj_exp
                    E_inj = 10**E_inj_exp
                case (13) ! shockless
                    if (arg == "true") then
                        shockless = .true.
                    else if (arg == "false") then
                        shockless = .false.
                    else
                        print *, "Invalid argument '", arg, "' for flag '", flag, "'."
                        call error("Invalid argument error", 0)
                    end if
                case (14) ! init-z-ax
                    if (arg == "true") then
                        init_z = .true.
                    else if (arg == "false") then
                        init_z = .false.
                    else
                        print *, "Invalid argument '", arg, "' for flag '", flag, "'."
                        call error("Invalid argument error", 0)
                    end if
                case (15) ! no_stepsize_corr
                    if (arg == "true") then
                        no_stepsize_corr = .true.
                    else if (arg == "false") then
                        no_stepsize_corr = .false.
                    else
                        print *, "Invalid argument", arg, " for flag '--shockless'."
                        call error("Invalid argument error", 0)
                    end if
                end select
                exit
            elseif (j == n_flags) then
                call error('Argument error, invalid flag: '//flag, 0)
            end if
        end do
    end do
end subroutine parse_cmd_args

subroutine init_general(myid, n_proc)
    use internal
    use user_variables
    use SNR_data

    implicit none
    integer, intent(in) :: myid, n_proc
    double precision :: v_EDST, v_shock

    ! Parse command line arguments and apply given settings/configuration
    call parse_cmd_args

    ! Add config metadata to the filename
    if (shockless) then
        filename = trim(basename_rw)
        if (no_stepsize_corr) then
            filename = filename//"_iso-stepsize"
        end if
        if (init_z) then
            filename = filename//"_init-z-ax"
        end if
        write (shockless_t_max_str, '(f10.3)') shockless_t_max
        filename = filename//"_t-max"//trim(adjustl(shockless_t_max_str))
    else
        if (inj_model == 0) then
            if (gamma_shock_set) then
                write (gamma_shock_str, '(f10.3)') gamma_shock
                filename = filename//"_gamma"//trim(adjustl(gamma_shock_str))
            else
                write (v_shock_str, '(f10.3)') v_shock(0)
                filename = filename//"_vshock"//trim(adjustl(v_shock_str))
            end if
        else if (inj_model == 1) then
            filename = filename//"_injmod1"
        else if (inj_model == 2) then
            filename = filename//"_injmod2"
        end if
    end if
    write (theta_max_str, '(f10.3)') pi*theta_max_pi_fraction
    write (n_sets_str, '(I10)') n_sets
    write (n_start_str, '(I10)') n_start
    write (n_proc_str, '(I10)') n_proc
    write (n_proc_str, '(I10)') n_proc
    write (E_inj_exp_str, '(f10.3)') log10(E_inj)
    filename = filename//'_theta-max'//trim(adjustl(theta_max_str))
    filename = filename//"_nsets"//trim(adjustl(n_sets_str))
    filename = filename//"_nstart"//trim(adjustl(n_start_str))
    filename = filename//"_E-inj-exp"//trim(adjustl(E_inj_exp_str))
    filename = filename//"_nproc"//trim(adjustl(n_proc_str))

    ! Print the filename with added metadata
    write (*, *) "Filename metadata: ", filename

    ! Allocate dynamic arrays

    !initialisation for random number (NumRec):
    iseed = 15321 + 2*(1 + iseed_shift)*(myid + 1)
    n_in = 0                         ! # of particles on stack
    n_out = 0                        ! # of escaped particles
    d_f = log10(E_min) - 0.1d0       ! init internal variables
    d_em = log10(E_min_em) - 0.1d0
    if (myid == 0) then
        !call test_int
        open (unit=99, file='Data/error'//filename)
        write (*, *) 'iseed = ', iseed
    end if

    !transition time between ED and ST phases (yr) (McKee)  ! init SNR_data variables
    t_ch = 423.d0/sqrt(E_snr/1.d51)*(M_ej/M_sun)**(5.d0/6.d0)/n_ISM**(1.d0/3.d0) ! yr
    R_ch = 3.07d0*(M_ej/M_sun)**(1.d0/3.d0)/n_ISM**(1.d0/3.d0)! pc
    v_ch = R_ch/t_ch

    t_EDST = 0.495d0*t_ch  ! yr
    R_EDST = 0.727d0*R_ch ! pc
    v_EDST = sqrt(2.d0*E_snr/M_ej)

    !start and end=transition ST-PDS, Truelove/McKee for zeta_m = 1.
    t_inj_fin = 1.33d4*(E_snr/1d51)**(3.d0/14.d0)/n_ISM**(4.d0/7.d0)
    !t_inj_fin = t_max
    if (myid == 0 .and. inj_model > 0) then
        write (*, *)
        write (*, *) 't_EDST = ', real(t_EDST), '/yr, v_EDST = ', real(v_EDST)
        write (*, *) 'v_ch = ', real(v_ch)
        write (*, *) 't_inj_fin = ', real(t_inj_fin), '/yr'
        write (*, *) 'compared to'
        write (*, *) 't_max     =', real(t_max), '/yr'
    end if
end subroutine init_general

subroutine init_inject_spec
    !acceptance-rejection method, comparing dN/dEdt with f(t)=K*t**alpha_f
    !with  f(t_EDST) = dN/dEdt(t_EDST) => K =  dN/dEdt(t_EDST)/t_EDST**alpha_f
    use SNR_data
    use internal; use result

    implicit none
    double precision :: dNdEdt

    t_inj_fin = t_max
    d_time_inj = log10(t_inj_fin/t_inj_init)/dble(n_tbin_in) ! linear or log??
    d_time_out = log10(t_max/t_inj_init)/dble(n_tbin_out)     ! linear or log??
    select case (inj_model)
    case (0)                     ! stationary
        alpha_f = 1.d0           ! does not matter
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
        x(3) = t_shock(t)
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
        !t = t_inj_init
        !t = 2.d3
        x(1) = 0.d0
        x(2) = 0.d0
        x(3) = t_shock(t)               ! delta function at shock, planar approx
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
