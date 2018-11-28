#!/bin/make

CC = icc
#OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__ 
OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__ -D__CHECKPOINT__ 
#OPTS = -g -fp-trap=common -mkl=sequential -D__SERIAL__ -D__FLOAT64__ -D__CHECKPOINT__ 
#OPTS = -O2 -Wall -mkl=sequential -D__SERIAL__ -D__FLOAT64__  
#OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__  
#CC = gcc

#OPTS = -O2 -D__SERIAL__ -D__FLOAT64__
SSALDIR = ../SSAL

NAME = ssal
#OPTS = -g   -D__SERIAL__ 
#OPTS = -g
EXPDIR=examples
SRCDIR=src

EXESRC = $(EXPDIR)/testfit_FKPP.c $(EXPDIR)/testfit_PF.c $(EXPDIR)/testfit_GPF.c $(EXPDIR)/ABC-Examples/ABCREJ_FKPP.c $(EXPDIR)/ABC-Examples/ABCREJ_FKPP_MULTDATA.c $(EXPDIR)/ABC-Examples/ABCREJ_PF_MULTDATA.c $(EXPDIR)/ABC-Examples/ABCREJ_GPF_MULTDATA.c $(EXPDIR)/ABC-Examples/ABCREJ_PF.c $(EXPDIR)/ABC-Examples/ABCREJ_GPF.c $(EXPDIR)/ABC-Examples/ABCREJ_VGPF.c $(EXPDIR)/ABC-Examples/ABCMCMC_GPF.c $(EXPDIR)/ABC-Examples/ABCMCMC_GPF_MULTDATA.c
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
