FC = mpifort
DFLAGS = -Wall -C -g -fbacktrace #-ffpe-trap=zero,invalid,overflow,underflow
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
			$(LIB_DIR)/acceleration.f90

# Main program -- development flags
test: $(LIB_FILES) $(SRC_DIR)/main.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/main.f90 -o $(BIN_DIR)/$@

# Main program -- production flags
main: $(LIB_FILES) $(SRC_DIR)/main.f90
	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/main.f90 -o $(BIN_DIR)/$@

# Generation of small angles (back rotation etc)
small_angle: $(LIB_FILES) $(SRC_DIR)/small_angle.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/small_angle.f90 -o $(BIN_DIR)/$@

# Print stuff for debugging etc. 
print_params: $(LIB_FILES) $(SRC_DIR)/print_parameters.f90
	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/print_parameters.f90 -o $(BIN_DIR)/$@


##############
# DEPRECATED #
##############

####################################################
## isotropic rw, L-trans advection and energy gain #
####################################################
#iso_test: $(LIB_FILES) $(SRC_DIR)/main_iso.f90
#	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/main_iso.f90 -o $(BIN_DIR)/$@
#
#iso: $(LIB_FILES) $(SRC_DIR)/main_iso.f90
#	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/main_iso.f90 -o $(BIN_DIR)/$@
#
#######################################################################
## isotropic rw, galilean advection and average energy gain per cycle #
#######################################################################
#old: $(LIB_FILES) $(SRC_DIR)/old_main.f90
#	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES) $(SRC_DIR)/old_main.f90 -o $(BIN_DIR)/$@
#
#old_test: $(LIB_FILES) $(SRC_DIR)/old_main.f90
#	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES) $(SRC_DIR)/old_main.f90 -o $(BIN_DIR)/$@
#
#########################################################
## Shock less random walks (isotropic and pitcha angle) #
#########################################################
#isotropic_rw: $(LIB_FILES_RW) $(SRC_DIR)/isotropic_rw.f90
#	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/isotropic_rw.f90 -o $(BIN_DIR)/$@
#
#small_angle_rw: $(LIB_FILES_RW) $(SRC_DIR)/small_angle_rw.f90
#	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/small_angle_rw.f90 -o $(BIN_DIR)/$@
#
#pitch_angle_rw: $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90
#	$(FC) -J $(MOD_DIR) $(PFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90 -o $(BIN_DIR)/$@
#
#pitch_angle_rw_test: $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90
#	$(FC) -J $(MOD_DIR) $(DFLAGS) $(LIB_FILES_RW) $(SRC_DIR)/pitch_angle_rw.f90 -o $(BIN_DIR)/$@
#
###################################
# Just generation of pitch angles #
###################################

#small_angle: $(SRC)/modules.f90 $(SRC)/small_angle.f90 $(SRC)/init.f90 $(SRC)/functions.f90 $(SRC)/output.f90 $(SRC)/aux.f90
#	gfortran -J$(MOD) -O -w $(SRC)/modules.f90 $(SRC)/small_angle.f90 $(SRC)/init.f90 $(SRC)/functions.f90 $(SRC)/output.f90 $(SRC)/aux.f90 -o $(BIN)/small_angle
#
#
