!File: output101.f90                      !

subroutine output(set, n_proc)
    use result; use SNR_data; use internal; use user_variables
    use constants, only: pi

    implicit none
    integer set, n_proc, i, j, pid
    double precision l, E, m, nu_tot, N_bin_j
    double precision :: bin_size, bin_size_aniso
    integer :: bin_num

    call banner(n_proc, set)
    n_tot = n_tot + n_start*n_proc
    write (*, *) 'set,n_tot', set, n_tot
    !open (20, file='Data/spec_aneut'//filename)    ! -8
    !open (21, file='Data/spec_aprot'//filename)    ! -7
    !open (22, file='Data/spec_anum'//filename)     ! -5
    !open (23, file='Data/spec_anue'//filename)     ! -4
    !open (24, file='Data/spec_pos'//filename)      ! -1
    !open (25, file='Data/spec_gam'//filename)      !  0
    !open (26, file='Data/spec_ele'//filename)
    !open (27, file='Data/spec_nue'//filename)
    !open (28, file='Data/spec_num'//filename)
    open (29, file=trim(outdir)//'/spec_prot'//filename)
    !open (30, file='Data/spec_neut'//filename)
    !open (50, file='Data/spec_nu'//filename)
    open (51, file=trim(outdir)//'/enumerate_spec_prot'//filename)
    open (60, file='Data/cross-angle-dist'//filename)
    open (70, file='Data/aniso-cross-angle-dist'//filename)
    do j = 1, n_enbin
        nu_tot = 0.d0
        pid = stable_pid(10) ! Proton id
        l = dble(j)*dn + d_f
        E = 10.d0**l
        m = En_f_tot(pid, j)/(dble(n_tot)*log(10.d0)*dn)
        if (m > 0.d0) write (19 + 10, 23) E, m, log10(E), log10(m)
        N_bin_j = exit_energy_enumerate_dist_tot(j)
        if (N_bin_j > 0.d0) then
            write (51, 23) E, N_bin_j, log10(E), log10(N_bin_j)
        end if
        ! With interactions on:
        !do i = 1, n_stable
        !    pid = stable_pid(i)
        !    l = dble(j)*dn + d_f
        !    E = 10.d0**l
        !    m = En_f_tot(pid, j)/(dble(n_tot)*log(10.d0)*dn)
        !    if (m > 0.d0) write (19 + i, 23) E, m, log10(E), log10(m)
        !    !if (m>0.d0) write(20+i,23) E,m,log10(E),log10(m)
        !    if (abs(pid) == 4 .or. abs(pid) == 5) nu_tot = nu_tot + m
        !end do
        !write (50, 23) E, nu_tot
    end do
    !close (20); close (21); close (22); close (23); close (24); close (25);
    !close (26); close (27); close (28); close (29); close (30); close (50);
    bin_size = pi/n_angle_bins
    do j = 1, n_angle_bins, 1
        write (60, 23) &
            j*bin_size, &
            cross_angle_distribution_first_tot(j), &
            cross_angle_distribution_updown_tot(j), &
            cross_angle_distribution_downup_tot(j)
    end do
    if (inj_model == 0 .and. gamma_shock_const >= 100) then
        bin_size_aniso = (pi/2.d0)/n_angle_bins_aniso
        do j = 1, n_angle_bins_aniso, 1
            write (70, 23) &
                j*bin_size_aniso, &
                cross_angle_distribution_aniso_updown_tot(j)
        end do
    end if
    close (29)
    close (51)
    close (60)
    close (70)

23  format(24E16.6)
    En_f_tot = 0.d0
    cross_angle_distribution_first_tot = 0.d0
    cross_angle_distribution_updown_tot = 0.d0
    cross_angle_distribution_downup_tot = 0.d0
end subroutine output

subroutine banner(n_proc, i)
    use user_variables
    use SNR_data
    use internal

    implicit none
    integer n_proc, i

    write (*, *)
    write (*, *) ' files saved as ', filename
    write (*, *)
    write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    if (i == 0) write (*, *) '!  SNR -- start of the main program      !'
    if (i == n_sets) write (*, *) '!  SNR - end of the main program         !'
    write (*, *) '!  # processes,           ', n_proc, '  !'
    write (*, *) '!  # MC sets,             ', n_sets, '  !'
    write (*, *) '!  injection model,       ', inj_model, '  !'
    !write(*,*) '!  fragmentation model,   ',fragm_model,'  !'
    write (*, *) '!  injected nuclei/set,   ', n_start*n_proc, '  !'
    if (xi_scat >= 1.d0) then
        write (*, *) '!  interaction enhanced by', int(xi_scat), '  !'
        write (*, *) '!  secondary spectra are rescaled        !'
    else
        call error('xi_scat<1 ..?', 1)
    end if
    write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end subroutine banner

subroutine init_mpi_io(comm_world, file_handle, fname, ierr)
    ! Opens mpi file for concurrent I/O
    use mpi_f08
    implicit none
    type(MPI_COMM), intent(in) :: comm_world
    type(MPI_FILE), intent(out) :: file_handle
    character(len=*), intent(in) :: fname
    integer, intent(out) :: ierr

    call MPI_FILE_OPEN( &
        MPI_COMM_WORLD, &
        fname, &
        MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
        MPI_INFO_NULL, &
        file_handle, &
        ierr &
        )

    if (ierr > 0) then
        call MPI_FILE_DELETE( &
            fname, &
            MPI_INFO_NULL, &
            ierr &
            )

        call MPI_FILE_OPEN( &
            MPI_COMM_WORLD, &
            fname, &
            MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_EXCL, &
            MPI_INFO_NULL, &
            file_handle, &
            ierr &
            )
    end if
end subroutine init_mpi_io
