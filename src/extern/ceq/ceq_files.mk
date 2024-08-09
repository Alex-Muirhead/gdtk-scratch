# object and source files needed to build ceq
# Note that we list just the D-language source files.
#
SRC_DIR := $(CEQ_DIR)/source
CEQ_SRC_FILES := $(SRC_DIR)/ceq.d

CEQ_OBJ_FILES := $(SRC_DIR)/ceq.o \
	$(SRC_DIR)/common.o \
	$(SRC_DIR)/linalg.o \
	$(SRC_DIR)/pt.o \
	$(SRC_DIR)/rhou.o \
	$(SRC_DIR)/ps.o \
	$(SRC_DIR)/rhot.o \
	$(SRC_DIR)/thermo.o
