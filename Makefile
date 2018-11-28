#!/bin/make

CC = icc
#OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__ 
OPTS = -g -mkl=sequential -D__SERIAL__ -D__FLOAT64__ -D__CHECKPOINT__ 
#OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__  
#CC = gcc

#OPTS = -O2 -D__SERIAL__ -D__FLOAT64__
SSALDIR = /home/n5685273/QUT/IF49/Projects/git_repos/SSAL

NAME = ssal
#OPTS = -g   -D__SERIAL__ 
#OPTS = -g
EXPDIR=examples
SRCDIR=src

EXESRC = $(EXPDIR)/ABC-Examples/MLABC_VarTest.c $(EXPDIR)/ABC-Examples/ABCvsMLMC_BCRN.c $(EXPDIR)/ABC-Examples/MLABC_TB.c $(EXPDIR)/ABC-Examples/ABCPCR_TB.c $(EXPDIR)/ABC-Examples/ABCMCMC_TB.c $(EXPDIR)/ABC-Examples/ABCREJ_TB.c $(EXPDIR)/ABC-Examples/ABC_toyProblem.c $(EXPDIR)/ABC-Examples/ABCMCMC_toyProblem.c $(EXPDIR)/ABC-Examples/ABCPCR_toyProblem.c $(EXPDIR)/ABC-Examples/MLABC_toyProblem.c $(EXPDIR)/ABC-Examples/ABCSMC_detLV.c $(EXPDIR)/ABC-Examples/ABCREJ_detLV.c $(EXPDIR)/ABC-Examples/MLABC_detLV.c $(EXPDIR)/ABC-Examples/ABCREJ_detGRN.c $(EXPDIR)/ABC-Examples/ABCSMC_stoLV.c $(EXPDIR)/ABC-Examples/MLABC_stoLV.c $(EXPDIR)/ABC-Examples/ABCREJ_stoLV.c $(EXPDIR)/MLMC-Examples/MLMC_Dimerisation.c
#SRC = $(SRCDIR)/sabccdfs.c $(SRCDIR)/smlabccdfs.c $(SRCDIR)/sabcrs.c $(SRCDIR)/sabcrjs.c $(SRCDIR)/dabccdfs.c $(SRCDIR)/dmlabccdfs.c $(SRCDIR)/dabcrs.c $(SRCDIR)/dabcrjs.c
SRC = $(SRCDIR)/mcl.c $(SRCDIR)/abc/dabccdfs.c $(SRCDIR)/abc/dmlabccdfs.c $(SRCDIR)/abc/dabcrs.c $(SRCDIR)/abc/dabcmcmc.c $(SRCDIR)/abc/dabcpcr.c $(SRCDIR)/abc/dabcrjs.c $(SRCDIR)/abc/dmlabcnls.c $(SRCDIR)/abc/dabcns.c $(SRCDIR)/mc/dmcint.c $(SRCDIR)/mc/dmcintd.c $(SRCDIR)/mc/dmcintv.c $(SRCDIR)/mc/dmlmcs.c $(SRCDIR)/mc/dmlmcnls.c $(SRCDIR)/data/copyDataset.c $(SRCDIR)/utils/monotone.c $(SRCDIR)/utils/cdf2pmf.c $(SRCDIR)/abc/dmlabcnlps.c $(SRCDIR)/pm/dmmh.c
OBJS = $(SRC:.c=.o)
EXEOBJS=$(EXESRC:.c=.o)
EXE = $(EXESRC:.c=)
INCDIR = $(SSALDIR)/include
LIBDIR = $(SSALDIR)/lib
INC = -I ./include -I $(INCDIR) 
LIBS = $(LIBDIR)/libssal.a -lm 

.SUFFIXES: .c .o

all: $(EXE) 
	@echo Binaries $(EXE) created!

.c.o: 
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

$(EXE): $(EXEOBJS) $(OBJS) 
	$(CC) $(OPTS) -o $@  $(INC) $@.o $(OBJS)   $(LIBS)
clean:
	set nonomatch; rm -f $(EXE) $(EXEOBJS) $(OBJS) 
