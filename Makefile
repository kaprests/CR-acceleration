SRC = ./src
MOD = ./Modules
BIN = ./bin

test: $(SRC)/modules101.f90 $(SRC)/main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90
	gfortran -J$(MOD) -Wall -C -g -fbacktrace $(SRC)/modules101.f90 $(SRC)/main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 -o $(BIN)/test


v101: modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90
	gfortran -JModules -O -w modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90 -o v101


old: modules101.f90 main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90
	gfortran -JModules -O -w modules101.f90 old_main101.f90 init101.f90 functions101.f90 output101.f90 aux101.f90 -o old


