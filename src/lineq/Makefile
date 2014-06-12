OBJECTS = gaussjbs.o test testinv

CC = gcc

TESTDEPS = ../misc/matrix.o ../lineq/gaussjbs.o
all: $(OBJECTS)

test: $(TESTDEPS)

testinv: $(TESTDEPS)

clean:
	rm $(OBJECTS)
