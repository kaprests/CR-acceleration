SRC = ./src
MOD = ./Modules
BIN = ./bin

test: $(SRC)/modules101.f90 $(SRC)/main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 $(SRC)/acceleration.f90
	gfortran -J$(MOD) -Wall -C -g -fbacktrace $(SRC)/modules101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 $(SRC)/acceleration.f90 $(SRC)/main101.f90 -o $(BIN)/test

v101: $(SRC)/modules101.f90 $(SRC)/main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90
	gfortran -J$(MOD) -O -w $(SRC)/modules101.f90 $(SRC)/main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 -o $(BIN)/v101

old: $(SRC)/modules101.f90 $(SRC)/old_main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90
	gfortran -J$(MOD) -O -w $(SRC)/modules101.f90 $(SRC)/old_main101.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 -o $(BIN)/old
