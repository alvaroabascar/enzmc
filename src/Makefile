BASE := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

CFLAGS = -Wall -g -I$(BASE)/include/

LDLIBS = -lm

CC = gcc

DEPS = montecarlo/montecarlo.o nlr/lvmrq.o random/random.o misc/mathlib.o misc/matrix.o lineq/gaussjbs.o

all:	enzmc

enzmc: $(DEPS)

objects:
	cd lineq; make
	cd misc; make
	cd models; make
	cd nlr; make
	cd random; make

clean:
	rm enzmc

cleanall:
	rm enzmc
	cd lineq; make clean 2>/dev/null
	cd misc; make clean 2>/dev/null
	cd models; make clean 2>/dev/null
	cd nlr; make clean 2>/dev/null
	cd random; make clean 2>/dev/null
	cd models; make clean 2>/dev/null
	cd montecarlo; make clean 2>/dev/null
