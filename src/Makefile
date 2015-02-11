CFLAGS = -Iinclude/ -g

LDLIBS = -lm

CC = gcc

DEPS = models/models.o montecarlo/montecarlo.o nlr/lvmrq.o random/random.o misc/mathlib.o misc/matrix.o lineq/gaussjbs.o

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
	cd lineq; make clean
	cd misc; make clean
	cd models; make clean
	cd nlr; make clean
	cd random; make clean
