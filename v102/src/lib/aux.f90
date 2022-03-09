! File: aux.f90
! Auxillary/utility subroutines and functions

subroutine error(string, s)
    !Error handling
    implicit none
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

function ran0()
    !random number generator from Numerical Recipes (Fortran90)
    use internal, only: iseed

    implicit none
    integer, parameter :: K4B = selected_int_kind(9)
    !integer(K4B), intent(inout) :: iseed
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
