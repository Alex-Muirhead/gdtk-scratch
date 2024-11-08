# makefile for kinetics module
# Can be used to perform unit tests
# and build stand-alone programs.

DMD := ldmd2

VPATH = kinetics

KINETICS_DIR := .

ifeq ($(shell uname -s), Darwin)
    PLATFORM := macosx
else
    PLATFORM := linux
endif
$(info PLATFORM=$(PLATFORM))

INSTALL_DIR ?= $(HOME)/gdtkinst
BUILD_DIR ?= ../../build

TEST_PROGRAMS := chemistry_update_test \
		equilibrium_update_test \
		rate_constant_test \
		reaction_test \
		reaction_mechanism_test \
		two_temperature_argon_kinetics_test \
		vib_specific_co_kinetics_test \
		vib_specific_co_mixture_kinetics_test \
		two_temperature_argon_with_ideal_gas_test \
		two_temperature_gasgiant_kinetics_test

NTYPES_DIR := ../ntypes
include $(NTYPES_DIR)/ntypes_files.mk

UTIL_DIR := ../util
include $(UTIL_DIR)/util_files.mk

NM_DIR := ../nm
include $(NM_DIR)/nm_files.mk

GAS_DIR := ../gas
include $(GAS_DIR)/gas_files.mk
LIBGASF := $(GAS_DIR)/libgasf.a

CEQ_DIR := ../extern/ceq
LIBCEQ := $(CEQ_DIR)/libceq.a
include $(CEQ_DIR)/ceq_files.mk

include kinetics_files.mk

LUA_DIR := ../../extern/lua-5.4.3
LIBLUA := ${LUA_DIR}/install/lib/liblua.a
LIBLUAPATH := ${LUA_DIR}/lib

DLINKFLAGS :=  -w -L-ldl
ifeq ($(PLATFORM), macosx)
    DLINKFLAGS += -L-ld_classic
endif

BIN_DIR := $(BUILD_DIR)/bin
LIB_DIR := $(BUILD_DIR)/lib

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

prep-chem: prep_chem.lua reaction.lua lex_elems.lua | $(BIN_DIR) $(LIB_DIR)
	cp $? $(LIB_DIR)
	cp $< $(BIN_DIR)/$@; chmod +x $(BIN_DIR)/$@

chemkin2eilmer: chemkin2eilmer.lua lex_elems.lua reaction.lua | $(BIN_DIR) $(LIB_DIR)
	cp $? $(LIB_DIR)
	cp $< $(BIN_DIR)/$@; chmod +x $(BIN_DIR)/$@

prep-kinetics: prep_kinetics.lua mechanism.lua lex_elems.lua | $(BIN_DIR) $(LIB_DIR)
	cp $? $(LIB_DIR)
	cp $< $(BIN_DIR)/$@; chmod +x $(BIN_DIR)/$@


test: $(TEST_PROGRAMS)
	- mv -n $(TEST_PROGRAMS) kinetics/.
	cd kinetics; tclsh kinetics-package-test.tcl

clean:
	- rm -f *.o *.mod *.obj
	- cd kinetics; rm -f $(TEST_PROGRAMS)
	- cd kinetics; rm -f two_temperature_argon_kinetics_test_results.data
	- rm -rf $(BUILD_DIR)/*
	- cd $(LUA_DIR); make clean
	- cd $(CEQ_DIR); make clean

$(LIBLUA):
	cd $(LUA_DIR); make $(PLATFORM) local PLATFORM=$(PLATFORM)

$(LIBCEQ):
	cd $(CEQ_DIR); make

$(TEST_PROGRAMS): \
	$(KINETICS_FILES) $(GAS_FILES) $(CEQ_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
	$(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
	$(LIBCEQ) $(LIBLUA)
	ldmd2 -of$@ -debug -g -dip1008 -version=$@ $^ $(DLINKFLAGS)
