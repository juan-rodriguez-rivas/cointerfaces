# Makefile.in for Easel's miniapps
# SVN $Id: Makefile.in 509 2010-02-07 22:56:55Z eddys $

top_srcdir = ../..
srcdir     = .
VPATH      =  ../..

SHELL       = /bin/sh
SUFFIXES    = 
SUFFIXES    = .c .o

## Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin


## your compiler and compiler flags.
#
CC        = gcc -std=gnu99
CFLAGS    = -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing
LDFLAGS   = 
SIMDFLAGS = -msse2
DEFS      = -DHAVE_CONFIG_H
LIBS      = -leasel   -lm

## where is Easel?
#
ESLDIR    = ..

## list of the miniapps to compile.
#
PROGS = esl-afetch\
	esl-alimanip\
	esl-alimap\
	esl-alimask\
	esl-alimerge\
	esl-alistat\
	esl-cluster\
	esl-compalign\
	esl-compstruct\
	esl-construct\
	esl-histplot\
	esl-mask\
	esl-reformat\
	esl-selectn\
	esl-seqstat\
	esl-sfetch\
	esl-shuffle\
	esl-ssdraw\
	esl-weight


all: 	
	cd ${ESLDIR}; ${MAKE} easel_lib
	${MAKE} easel_progs

dev:
	cd ${ESLDIR}; ${MAKE} easel_lib
	${MAKE} easel_progs

check:
	cd ${ESLDIR}; ${MAKE} easel_lib
	${MAKE} easel_progs

easel_progs: ${PROGS}

${PROGS}: %: %.o  ${ESLDIR}/libeasel.a
	${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} -L${ESLDIR} -o $@ $@.o ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} -I. -I${ESLDIR} -I${srcdir} -I${srcdir}/.. -o $@ -c $<

install:
	cp ${PROGS} ${DESTDIR}${bindir}/

uninstall:
	for file in $(PROGS); do\
	   rm ${DESTDIR}${bindir}/$$file;\
	done

clean:
	-rm -f ${PROGS}
	-rm -f *.o *~
	-rm -f TAGS

distclean:
	${MAKE} clean
	-rm -f ${PROGS}
	-rm Makefile

TAGS:
	etags *.[ch] Makefile.in


################################################################
# Easel - a library of C functions for biological sequence analysis
# Version h3.0; March 2010
# Copyright (C) 2010 Howard Hughes Medical Institute.
# Other copyrights also apply. See the COPYRIGHT file for a full list.
# 
# Easel is distributed under the Janelia Farm Software License, a BSD
# license. See the LICENSE file for more details.
################################################################
