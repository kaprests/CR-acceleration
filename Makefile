FC = gfortran
DFLAGS = -Wall -C -g -fbacktrace
PFLAGS = -O -w

SRC_DIR = ./src
LIB_DIR = ./src/lib
MOD_DIR = ./Modules
BIN_DIR = ./bin

#LIB_FILES = $(wildcard $(LIB_DIR)/*.f90)
LIB_FILES = $(LIB_DIR)/modules101.f90 $(LIB_DIR)/init101.f90 $(LIB_DIR)/functions101.f90 $(LIB_DIR)/output101.f90 $(LIB_DIR)/aux101.f90 $(LIB_DIR)/acceleration.f90


test: $(LIB_FILES) $(SRC_DIR)/main101.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/main101.f90 -o $(BIN_DIR)/$@


old: $(LIB_FILES) $(SRC_DIR)/old_main101.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/old_main101.f90 -o $(BIN_DIR)/$@


main: $(LIB_FILES) $(SRC_DIR)/main101.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/main101.f90 -o $(BIN_DIR)/$@


small_angle: $(LIB_FILES) $(SRC_DIR)/small_angle.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/small_angle.f90 -o $(BIN_DIR)/$@


isotropic_rw: $(LIB_FILES) $(SRC_DIR)/isotropic_rw.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/isotropic_rw.f90 -o $(BIN_DIR)/$@

small_angle_rw: $(LIB_FILES) $(SRC_DIR)/small_angle_rw.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/small_angle_rw.f90 -o $(BIN_DIR)/$@

pitch_angle_rw: $(LIB_FILES) $(SRC_DIR)/pitch_angle_rw.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/pitch_angle_rw.f90 -o $(BIN_DIR)/$@

#small_angle: $(SRC)/modules101.f90 $(SRC)/small_angle.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90
#	gfortran -J$(MOD) -O -w $(SRC)/modules101.f90 $(SRC)/small_angle.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 -o $(BIN)/small_angle
#
#
