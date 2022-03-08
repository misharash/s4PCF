## MAKEFILE FOR gaussian_s4PCF. This compiles the gs4PCF.cpp file into the ./gs4PCF exececutable.

CC = gcc
CFLAGS = -g -O3 -Wall -MMD
CXXFLAGS = -DOPENMP -O3 -Wall -MMD

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math

EXEC	= gs4PCF
OBJS	= gs4PCF.o ./cubature/hcubature.o ./ransampl/ransampl.o
DEPS	= ${OBJS:.o=.d}

LD	= g++
LFLAGS	= -lgsl -lgslcblas -lgomp

.PHONY: main clean

main: $(EXEC)

$(EXEC):	$(OBJS) Makefile
	$(LD) $(OBJS) $(LFLAGS) -o $(EXEC)

clean:
	rm -f $(EXEC) $(OBJS) ${DEPS}

$(OBJS): Makefile
-include ${DEPS}
