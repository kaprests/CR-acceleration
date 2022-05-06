FC = mpifort
FC_NOMPI = gfortran
DFLAGS = -Wall -C -g -fbacktrace
PFLAGS = -O -w

SRC_DIR = ./src
LIB_DIR = ./src/lib
MOD_DIR = ./Modules
BIN_DIR = ./bin

LIB_FILES = $(LIB_DIR)/modules.f90\
			$(LIB_DIR)/init.f90\
			$(LIB_DIR)/functions.f90\
			$(LIB_DIR)/output.f90\
			$(LIB_DIR)/aux.f90\
			$(LIB_DIR)/random_walk.f90

LIB_FILES_TEST = $(LIB_DIR)/modules.f90\
			$(LIB_DIR)/init.f90\
			$(LIB_DIR)/functions.f90\
			$(LIB_DIR)/aux.f90\
			$(LIB_DIR)/random_walk.f90

all: dev main test

dev: $(LIB_FILES) $(SRC_DIR)/main.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/main.f90 -o $(BIN_DIR)/$@

main: $(LIB_FILES) $(SRC_DIR)/main.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/main.f90 -o $(BIN_DIR)/$@

test: $(LIB_FILES_TEST) $(SRC_DIR)/test.f90
	$(FC_NOMPI) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES_TEST) $(SRC_DIR)/test.f90 -o $(BIN_DIR)/$@

clean:
	rm -f fort.99 ./bin/main ./bin/dev ./bin/test
