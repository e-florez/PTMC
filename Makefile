# Makefile.in - Template for Fortran 2003 code Makefile

# Variables
SRC_DIR := src
MAIN := $(SRC_DIR)/atomic_PTMC_look-up_Table_v3.f03
SUBS := $(SRC_DIR)/DrawMolecularConfiguration_v9.f03 $(SRC_DIR)/ClockSeed_v1.f03 $(SRC_DIR)/RandomNumber_v1.f03 $(SRC_DIR)/DisplaceAtoms_v1.f03 $(SRC_DIR)/RotateCluster_v1.f03 $(SRC_DIR)/AtomicMass_v1.f03
OUT := PTMC.exe
MULTIHIST_SRC := $(SRC_DIR)/multihist_dCv.f
MULTIHIST_OUT := multihist_dCv.exe

# Compiler
FC := gfortran

# Default flags
FFLAGS := -O3 -fimplicit-none

# Check the command line argument and set the appropriate flags
ifeq ($(MAKECMDGOALS),openmp)
	FFLAGS += -fopenmp
endif
ifeq ($(MAKECMDGOALS),array)
	FFLAGS += -fbounds-check -fbacktrace -fcheck=all
endif
ifeq ($(MAKECMDGOALS),debug)
	FFLAGS += -fbounds-check -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Og
endif
ifeq ($(MAKECMDGOALS),debug-openmp)
	FFLAGS += -fopenmp -fbounds-check -g -Wall -Wextra -Warray-temporaries -Wconversion -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan -Og
endif

.PHONY: all openmp array debug debug-openmp clean

all: $(OUT) $(MULTIHIST_OUT)

openmp: $(OUT) $(MULTIHIST_OUT)

array: $(OUT) $(MULTIHIST_OUT)

debug: $(OUT) $(MULTIHIST_OUT)

debug-openmp: $(OUT) $(MULTIHIST_OUT)

$(OUT): $(MAIN) $(SUBS)
	$(FC) $(FFLAGS) $^ -o $@

$(MULTIHIST_OUT): $(MULTIHIST_SRC)
	$(FC) $< -ffpe-summary=none -o $@

clean:
	rm -f $(OUT) $(MULTIHIST_OUT)

