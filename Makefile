FC = mpifort
DFLAGS = -Wall -C -g -fbacktrace #-ffpe-trap=zero,invalid,overflow,underflow
PFLAGS = -O -w

SRC_DIR = ./src
LIB_DIR = ./src/lib
MOD_DIR = ./Modules
BIN_DIR = ./bin

#LIB_FILES = $(wildcard $(LIB_DIR)/*.f90)
LIB_FILES = $(LIB_DIR)/modules101.f90 $(LIB_DIR)/init101.f90 $(LIB_DIR)/functions101.f90 $(LIB_DIR)/output101.f90 $(LIB_DIR)/aux101.f90 $(LIB_DIR)/acceleration.f90

LIB_FILES_RW = $(LIB_DIR)/modules101.f90 $(LIB_DIR)/init101.f90 $(LIB_DIR)/functions101.f90 $(LIB_DIR)/output101.f90 $(LIB_DIR)/aux101.f90 

########################################################
# pitch angle, L-transformed advection and energy gain #
########################################################
test: $(LIB_FILES) $(SRC_DIR)/main.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/main.f90 -o $(BIN_DIR)/$@

main: $(LIB_FILES) $(SRC_DIR)/main.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/main.f90 -o $(BIN_DIR)/$@

###################################################
# isotropic rw, L-trans advection and energy gain #
###################################################
iso_test: $(LIB_FILES) $(SRC_DIR)/main_iso.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/main_iso.f90 -o $(BIN_DIR)/$@

iso: $(LIB_FILES) $(SRC_DIR)/main_iso.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/main_iso.f90 -o $(BIN_DIR)/$@

######################################################################
# isotropic rw, galilean advection and average energy gain per cycle #
######################################################################
old: $(LIB_FILES) $(SRC_DIR)/old_main101.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/old_main101.f90 -o $(BIN_DIR)/$@

old_test: $(LIB_FILES) $(SRC_DIR)/old_main101.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/old_main101.f90 -o $(BIN_DIR)/$@

########################################################
# Shock less random walks (isotropic and pitcha angle) #
########################################################
isotropic_rw: $(LIB_FILES_RW) $(SRC_DIR)/isotropic_rw.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/isotropic_rw.f90 -o $(BIN_DIR)/$@

small_angle_rw: $(LIB_FILES_RW) $(SRC_DIR)/small_angle_rw.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/small_angle_rw.f90 -o $(BIN_DIR)/$@

pitch_angle_rw: $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90 -o $(BIN_DIR)/$@

pitch_angle_rw_test: $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90 -o $(BIN_DIR)/$@

###################################
# Just generation of pitch angles #
###################################
small_angle: $(LIB_FILES) $(SRC_DIR)/small_angle.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/small_angle.f90 -o $(BIN_DIR)/$@

#small_angle: $(SRC)/modules101.f90 $(SRC)/small_angle.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90
#	gfortran -J$(MOD) -O -w $(SRC)/modules101.f90 $(SRC)/small_angle.f90 $(SRC)/init101.f90 $(SRC)/functions101.f90 $(SRC)/output101.f90 $(SRC)/aux101.f90 -o $(BIN)/small_angle
#
#
