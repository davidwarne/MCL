#!/bin/make

CC = gcc
#OPTS = -O2 -std=gnu99 -D__SERIAL__ -D__FLOAT64__ 
OPTS = -pg -g -std=gnu99 -D__SERIAL__ -D__FLOAT64__

NAME = ssal
#OPTS = -g   -D__SERIAL__ 
#OPTS = -g
EXPDIR=examples
SRCDIR=src

EXESRC = $(EXPDIR)/testMLABC.c
#SRC = $(SRCDIR)/sabccdfs.c $(SRCDIR)/smlabccdfs.c $(SRCDIR)/sabcrs.c $(SRCDIR)/sabcrjs.c $(SRCDIR)/dabccdfs.c $(SRCDIR)/dmlabccdfs.c $(SRCDIR)/dabcrs.c $(SRCDIR)/dabcrjs.c
SRC = $(SRCDIR)/dabccdfs.c $(SRCDIR)/dmlabccdfs.c $(SRCDIR)/dabcrs.c $(SRCDIR)/dabcrjs.c $(SRCDIR)/dmlabcnls.c $(SRCDIR)/dabcns.c $(SRCDIR)/dmcint.c $(SRCDIR)/dmcintd.c $(SRCDIR)/copyDataset.c
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