# Makefile for DaliLite.v5
# # Copyright (C) 2018 Liisa Holm.
# 
# ##############################################################################################
# # start of configuration section                                                             #
# ##############################################################################################
 SHELL = /bin/sh
 OPENMPI_PATH = /usr/lib64/openmpi/bin
 SRC = ../src
 BIN = ./
 CMBI = $(SRC)/dsspcmbi
# Compilers
 CC = gcc
 FC = gfortran
 CC_MPI = $(OPENMPI_PATH)/mpicc                   # C compiler
 FC_MPI = $(OPENMPI_PATH)/mpif90              # Fortran-90 compiler
 OPT = -O3                  # uncomment this for optimization
# Compiler options
 DEB = # -g                 # uncomment this for debugging
 DEFS =                     # place other -D types of things here
 CFLAGS = $(OPT) $(DEB) $(DEFS)
 LFLAGS =
# ##############################################################################################
# # end of configuration section                                                               #
# ##############################################################################################
# # objects
 OBJS = AccSurf.o CalcAccSurf.o  Contacts.o  Date.o  DsspCMBI.o  p2clib.o Vector.o 
#
 EXES = $(BIN)/serialcompare $(BIN)/pipe $(BIN)/filter95fitz $(BIN)/selfdccp $(BIN)/pipedccp $(BIN)/fssp $(BIN)/puu $(BIN)/puutos $(BIN)/dsspcmbi
#
 MODS =  dalicon.mod  dp.mod  filter95.mod  parsi.mod  soap.mod  wolf.mod
#
# # Top-level targets
 all: 
	rm -rf AccSurf.o CalcAccSurf.o  Contacts.o  Date.o  DsspCMBI.o  p2clib.o Vector.o  .//serialcompare .//pipe .//filter95fitz .//selfdccp .//pipedccp .//fssp .//puu .//puutos .//dsspcmbi .//mpicompare dalicon.mod  dp.mod  filter95.mod  parsi.mod  soap.mod  wolf.mod

	cc -c  -O ../src/dsspcmbi/AccSurf.c
	cc -c  -O ../src/dsspcmbi/CalcAccSurf.c
	cc -c  -O ../src/dsspcmbi/Contacts.c
	cc -c  -O ../src/dsspcmbi/Date.c
	cc -c  -O ../src/dsspcmbi/DsspCMBI.c
	cc -c  -O ../src/dsspcmbi/p2clib.c
	cc -c  -O ../src/dsspcmbi/Vector.c
	cc -o .//dsspcmbi DsspCMBI.o p2clib.o Date.o Vector.o AccSurf.o CalcAccSurf.o Contacts.o -lm
	gfortran -O3                     ../src/util.f ../src/subfitzfast.f ../src/parsi-admin.f ../src/parsi-align.f ../src/parsi-score.f ../src/parsi-stack.f ../src/parsi-1.f ../src/gagadccp.f ../src/gagagene.f ../src/gagahash.f ../src/gagatool.f ../src/lean.f ../src/clean.f ../src/ssap.f ../src/dimple.f ../src/triplet.f ../src/testi.f ../src/soap4.f ../src/wolf_original.f ../src/u3b-8.f  ../src/comparemodules.f ../src/serialcompare.f -o .//serialcompare
	gfortran -O3                     ../src/pipe.f ../src/util.f -o .//pipe
	gfortran -O3                     ../src/filter95fitz.f ../src/util.f ../src/u3b-8.f ../src/subfitzfast.f -o .//filter95fitz
	gfortran -O3                     ../src/selfdccp.f ../src/util.f -o .//selfdccp
	gfortran -O3                     ../src/pipedccp.f -o .//pipedccp
	gfortran -O3                     ../src/fssp.f ../src/util.f ../src/u3b-8.f -o .//fssp
	gfortran -O3                     ../src/main.f ../src/puu.f ../src/util.f -o .//puu
	gfortran -O3                     ../src/puutos.f ../src/util.f -o .//puutos
	$(FC) $(CFLAGS) $(SRC)/pipe.f $(SRC)/util.f -o $(BIN)/pipe
	$(FC) $(CFLAGS) $(SRC)/filter95fitz.f $(SRC)/util.f $(SRC)/u3b-8.f $(SRC)/subfitzfast.f -o $(BIN)/filter95fitz
	$(FC) $(CFLAGS) $(SRC)/selfdccp.f $(SRC)/util.f -o $(BIN)/selfdccp
	$(FC) $(CFLAGS) $(SRC)/pipedccp.f -o $(BIN)/pipedccp
	$(FC) $(CFLAGS) $(SRC)/fssp.f $(SRC)/util.f $(SRC)/u3b-8.f -o $(BIN)/fssp
	$(FC) $(CFLAGS) $(SRC)/main.f $(SRC)/puu.f $(SRC)/util.f -o $(BIN)/puu
	$(FC) $(CFLAGS) $(SRC)/puutos.f $(SRC)/util.f -o $(BIN)/puutos

 parallel: all $(BIN)/mpicompare
 
 install: clean all 

# programs
$(BIN)/mpicompare: 
	$(FC_MPI) $(CFLAGS) $(SRC)/util.f $(SRC)/subfitzfast.f $(SRC)/parsi-admin.f $(SRC)/parsi-align.f $(SRC)/parsi-score.f $(SRC)/parsi-stack.f $(SRC)/parsi-1.f $(SRC)/gagadccp.f $(SRC)/gagagene.f $(SRC)/gagahash.f $(SRC)/gagatool.f $(SRC)/lean.f $(SRC)/clean.f $(SRC)/ssap.f $(SRC)/dimple.f $(SRC)/triplet.f $(SRC)/testi.f $(SRC)/soap4.f $(SRC)/wolf_original.f $(SRC)/u3b-8.f  $(SRC)/comparemodules.f $(SRC)/mpicompare.f -o $(BIN)/mpicompare

$(BIN)/serialcompare:  
	$(FC) $(CFLAGS) $(SRC)/util.f $(SRC)/subfitzfast.f $(SRC)/parsi-admin.f $(SRC)/parsi-align.f $(SRC)/parsi-score.f $(SRC)/parsi-stack.f $(SRC)/parsi-1.f $(SRC)/gagadccp.f $(SRC)/gagagene.f $(SRC)/gagahash.f $(SRC)/gagatool.f $(SRC)/lean.f $(SRC)/clean.f $(SRC)/ssap.f $(SRC)/dimple.f $(SRC)/triplet.f $(SRC)/testi.f $(SRC)/soap4.f $(SRC)/wolf_original.f $(SRC)/u3b-8.f  $(SRC)/comparemodules.f $(SRC)/serialcompare.f -o $(BIN)/serialcompare

$(BIN)/pipe: 
	$(FC) $(CFLAGS) $(SRC)/pipe.f $(SRC)/util.f -o $(BIN)/pipe

$(BIN)/filter95fitz: 
	$(FC) $(CFLAGS) $(SRC)/filter95fitz.f $(SRC)/util.f $(SRC)/u3b-8.f $(SRC)/subfitzfast.f -o $(BIN)/filter95fitz

$(BIN)/selfdccp: 
	$(FC) $(CFLAGS) $(SRC)/selfdccp.f $(SRC)/util.f -o $(BIN)/selfdccp

$(BIN)/pipedccp: 
	$(FC) $(CFLAGS) $(SRC)/pipedccp.f -o $(BIN)/pipedccp

$(BIN)/fssp: 
	$(FC) $(CFLAGS) $(SRC)/fssp.f $(SRC)/util.f $(SRC)/u3b-8.f -o $(BIN)/fssp

$(BIN)/puu: 
	$(FC) $(CFLAGS) $(SRC)/main.f $(SRC)/puu.f $(SRC)/util.f -o $(BIN)/puu

$(BIN)/puutos: 
	$(FC) $(CFLAGS) $(SRC)/puutos.f $(SRC)/util.f -o $(BIN)/puutos

# dsspcmbi

DsspCMBI.o: 
	cc -c  -O $(CMBI)/DsspCMBI.c

p2clib.o:
	cc -c  -O  $(CMBI)/p2clib.c

Date.o:
	cc -c  -O  $(CMBI)/Date.c

Vector.o:
	cc -c  -O  $(CMBI)/Vector.c

AccSurf.o:
	cc -c  -O  $(CMBI)/AccSurf.c

CalcAccSurf.o:
	cc -c  -O  $(CMBI)/CalcAccSurf.c

Contacts.o:
	cc -c  -O  $(CMBI)/Contacts.c

$(BIN)/dsspcmbi: DsspCMBI.o p2clib.o Date.o Vector.o AccSurf.o CalcAccSurf.o Contacts.o
	cc -o $(BIN)/dsspcmbi DsspCMBI.o p2clib.o Date.o Vector.o AccSurf.o CalcAccSurf.o Contacts.o -lm


# clean slate
clean: 
	rm -rf $(OBJS) $(EXES) $(BIN)/mpicompare $(MODS)

 
