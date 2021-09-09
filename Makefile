test: modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90
	gfortran -JModules -Wall -C -g -fbacktrace modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90 -o test


v101: modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90
	gfortran -JModules -O -w modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90 -o v101


dev: modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90
	gfortran -JModules -O -w modules101.f90 devmain.f90 init101.f90 functions101.f90 output101.f90 aux101.f90 -o dev


