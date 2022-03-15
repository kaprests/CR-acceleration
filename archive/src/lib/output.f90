!=============================================================================!
!            file ~/Programme/SNR/Nuclei11/output101.f90                      !
!=============================================================================!
subroutine output(set, n_proc)
    use result; use SNR_data; use internal; use user_variables
    implicit none
    integer set, n_proc, j, pid
    double precision l, E, m, nu_tot, N_ltE_esc

    call banner(n_proc, set)

    n_tot = n_tot + n_start*n_proc
    write (*, *) 'set,n_tot', set, n_tot

    !open (20, file=trim(outdir)//'/spec_aneut'//filename)   ! -8
    !open (21, file=trim(outdir)//'/spec_aprot'//filename)   ! -7
    !open (22, file=trim(outdir)//'/spec_anum'//filename)    ! -5
    !open (23, file=trim(outdir)//'/spec_anue'//filename)    ! -4
    !open (24, file=trim(outdir)//'/spec_pos'//filename)     ! -1
    !open (25, file=trim(outdir)//'/spec_gam'//filename)     !  0
    !open (26, file=trim(outdir)//'/spec_ele'//filename)
    !open (27, file=trim(outdir)//'/spec_nue'//filename)
    !open (28, file=trim(outdir)//'/spec_num'//filename)
    open (29, file=trim(outdir)//'/spec_prot'//filename)
    !open (30, file=trim(outdir)//'/spec_neut'//filename)
    !open (50, file=trim(outdir)//'/spec_nu'//filename)
    open (60, file=trim(outdir)//'/esc_prot'//filename)      ! sensible unit number?

    do j = 1, n_enbin
        nu_tot = 0.d0
        l = dble(j)*dn + d_f
        E = 10.d0**l
        !do i = 1, n_stable
        !   pid = stable_pid(i)
        !   m = En_f_tot(pid, j)/(dble(n_tot)*log(10.d0)*dn)
        !   if (m > 0.d0) write (19 + i, 23) E, m, log10(E), log10(m)
        !   !if (m>0.d0) write(20+i,23) E,m,log10(E),log10(m)
        !   if (abs(pid) == 4 .or. abs(pid) == 5) nu_tot = nu_tot + m
        !end do

        ! Write proton energies/spectrum
        pid = stable_pid(10) ! proton
        m = En_f_tot(pid, j)/(dble(n_tot)*log(10.d0)*dn)
        if (m > 0.d0) write (19 + 10, 23) E, m, log10(E), log10(m)
        if (abs(pid) == 4 .or. abs(pid) == 5) nu_tot = nu_tot + m

        ! Write escape probability
        if (j > 1) then
            N_ltE_esc = sum(NE_esc_tot(1:j - 1))
        else
            N_ltE_esc = 0
        end if
        ! E, log10(E), p_esc(estimate)
        if (n_tot - N_ltE_esc > 0) write (60, 23) E, log10(E), NE_esc(j)/(n_tot - N_ltE_esc)
        !write (50, 23) E, nu_tot
    end do

    !close (20); close (21); close (22); close (23); close (24); close (25);
    !close (26); close (27); close (28); close (29); close (30); close (50);
    close (29)
    close (60); 
23  format(24E16.6)

    En_f_tot = 0.d0
end subroutine output
!============================================================================!
!============================================================================!
!subroutine output_finish
!   use user_variables, only: filename, outdir
!   use result
!   implicit none
!     open (9, file=trim(outdir)//'/num_crossings'//filename, form='unformatted')
!     write (9) num_crossings_total
!     close (9)
!
!     open (9, file=trim(outdir)//'/trajectories'//filename, form='unformatted')
!     write (9) trajectories
!     close (9)
!
!   deallocate(drift_distances)
!   deallocate(final_positions)
!
!   deallocate(exit_energies)
!   deallocate(num_crossings_total)
!   deallocate(trajectories)
!   deallocate(sample_positions)
!
!   open (10, file=trim(outdir_raw)//'/exit_energies'//filename//'.dat', form='unformatted')
!   write (10) exit_energies
!   close (10)
!end subroutine output_finish
!============================================================================!
!============================================================================!
subroutine banner(n_proc, i)
    use user_variables
    use SNR_data
    use internal
    use result, only: rel_energy_gain_total_sum
    implicit none
    integer n_proc, i
    double precision :: rel_energy_gain_total_average

    rel_energy_gain_total_average = rel_energy_gain_total_sum/(n_sets*n_start)

    write (*, *)
    write (*, *) ' files saved as ', filename
    write (*, *)
    write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    if (i == 0) write (*, *) '!  SNR -- start of the main program      !'
    if (i == n_sets) write (*, *) '!  SNR - end of the main program         !'
    write (*, *) '!  # processes,           ', n_proc, '  !'
    write (*, *) '!  # MC sets,             ', n_sets, '  !'
    write (*, *) '!  injection model,       ', inj_model, '  !'
!  write(*,*) '!  fragmentation model,   ',fragm_model,'  !'
    write (*, *) '!  injected nuclei/set,   ', n_start*n_proc, '  !'
    if (xi_scat >= 1.d0) then
        write (*, *) '!  interaction enhanced by', int(xi_scat), '  !'
        write (*, *) '!  secondary spectra are rescaled        !'
    else
        call error('xi_scat<1 ..?', 1)
    end if
!   write (*, *) '!  Average relative E gain, ', rel_energy_gain_total_average, '  !'
!   write (*, *) '!  Vshock,              ', shock_velocity, '  !'
!   write (*, *) '!  Vshock^2,            ', shock_velocity**2, '  !'
    write (*, *) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

end subroutine banner

!!!!!!!!!!!!
!! MPI IO !!
!!!!!!!!!!!!

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

!! Does not work -- can't pass allocatable array buffer
!subroutine write_mpi_io(file_handle, data_buf, myid, set, ierr)
!   ! Writes to mpi file (from given proc). Writes double precision floats.
!   use mpi_f08
!   use user_variables, only: n_sets
!   use result
!
!   implicit none
!   type(MPI_FILE), intent(in) :: file_handle
!   double precision, dimension(..), intent(out) :: data_buf
!   integer, intent(in) :: myid, set
!   integer, intent(out) :: ierr
!
!   integer :: data_buf_count, data_buf_byte_size
!   integer(kind=MPI_OFFSET_KIND) :: offset
!
!   data_buf_count = size(data_buf)
!   data_buf_byte_size = sizeof(data_buf)
!   offset = myid * n_sets * data_buf_byte_size + data_buf_byte_size * (set - 1)
!
!   print *, "!!!!!!!!!!!!!!!!"
!   print *, "subroutine data: "
!   print *, data_buf_count
!   print *, data_buf_byte_size
!   !print *, data_buf
!   print *, "!!!!!!!!!!!!!!!!"
!
!   print *, "WRITING DATA NAOOOOO"
!   call MPI_FILE_WRITE_AT(&
!      file_handle, &
!      offset, &
!      data_buf, &
!      data_buf_count, &
!      MPI_DOUBLE_PRECISION, &
!      MPI_STATUS_IGNORE, &
!      ierr &
!   )
!end subroutine write_mpi_io
