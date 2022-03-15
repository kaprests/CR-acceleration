program acceleration
   use mpi_f08

   implicit none
   integer myid, n_proc, ierr

   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,ierr)
   
   write(*, *) "Printing from proc: ", myid
   call MPI_FINALIZE(ierr)
end program acceleration
