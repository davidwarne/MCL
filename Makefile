#!/bin/make

CC = icc
#OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__ 
OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__ -D__CHECKPOINT__ 
#OPTS = -O2 -mkl=sequential -D__SERIAL__ -D__FLOAT64__  
#CC = gcc
#OPTS = -O2 -D__SERIAL__ -D__FLOAT64__

NAME = ssal
#OPTS = -g   -D__SERIAL__ 
#OPTS = -g
EXPDIR=examples
SRCDIR=src

EXESRC = $(EXPDIR)/MLABC_VarTest.c $(EXPDIR)/ABCvsMLMC_BCRN.c $(EXPDIR)/MLABC_TB.c $(EXPDIR)/ABCPCR_TB.c $(EXPDIR)/ABCMCMC_TB.c $(EXPDIR)/ABCREJ_TB.c $(EXPDIR)/ABC_toyProblem.c $(EXPDIR)/ABCMCMC_toyProblem.c $(EXPDIR)/ABCPCR_toyProblem.c $(EXPDIR)/MLABC_toyProblem.c $(EXPDIR)/ABCSMC_detLV.c $(EXPDIR)/ABCREJ_detLV.c $(EXPDIR)/MLABC_detLV.c $(EXPDIR)/ABCREJ_detGRN.c $(EXPDIR)/ABCSMC_stoLV.c $(EXPDIR)/MLABC_stoLV.c $(EXPDIR)/ABCREJ_stoLV.c
#SRC = $(SRCDIR)/sabccdfs.c $(SRCDIR)/smlabccdfs.c $(SRCDIR)/sabcrs.c $(SRCDIR)/sabcrjs.c $(SRCDIR)/dabccdfs.c $(SRCDIR)/dmlabccdfs.c $(SRCDIR)/dabcrs.c $(SRCDIR)/dabcrjs.c
SRC = $(SRCDIR)/libabc.c $(SRCDIR)/dabccdfs.c $(SRCDIR)/dmlabccdfs.c $(SRCDIR)/dabcrs.c $(SRCDIR)/dabcmcmc.c $(SRCDIR)/dabcpcr.c $(SRCDIR)/dabcrjs.c $(SRCDIR)/dmlabcnls.c $(SRCDIR)/dabcns.c $(SRCDIR)/dmcint.c $(SRCDIR)/dmcintd.c $(SRCDIR)/copyDataset.c $(SRCDIR)/monotone.c $(SRCDIR)/cdf2pmf.c $(SRCDIR)/dmlabcnlps.c
OBJS = $(SRC:.c=.o)
EXEOBJS=$(EXESRC:.c=.o)
EXE = $(EXESRC:.c=)
INCDIR = /home/n5685273/QUT/IF49/Projects/BaselineSSA/SSAL/include
LIBDIR = /home/n5685273/QUT/IF49/Projects/BaselineSSA/SSAL/lib
INC = -I ./include -I $(INCDIR) 
LIBS = $(LIBDIR)/libssal.a -lm 

.SUFFIXES: .c .o

all: $(EXE) 
	@echo Binaries $(EXE) created!

.c.o: 
	$(CC) $(OPTS) $(PROFILE) -c $< -o $@ $(INC) 

#$(BIN): $(OBJS)
#	$(CC) $(OPTS) $(PROFILE)  $(OBJS) -o $(BIN) -lm 
#	@echo Binary created!!

$(EXE): $(EXEOBJS) $(OBJS) 
	$(CC) $(OPTS) -o $@  $(INC) $@.o $(OBJS)   $(LIBS)
clean:
	set nonomatch; rm -f $(EXE) $(EXEOBJS) $(OBJS) 
