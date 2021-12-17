CC = gcc
CFLAGS = -g -O3 -Wall -MMD
CXXFLAGS = -DOPENMP -O3 -Wall -MMD

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math

EXEC	= s4PCF
SRC	= $(wildcard *.cpp)
OBJS	= ${SRC:.cpp=.o}
DEPS	= ${OBJS:.o=.d}

LD	= g++
#LFLAGS	= -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lgomp
LFLAGS	= -lgomp

.PHONY: main clean

main: $(EXEC)

$(EXEC):	$(OBJS) Makefile
	$(LD) $(OBJS) $(LFLAGS) -o $(EXEC)

clean:
	rm -f $(EXEC) $(OBJS) ${DEPS}

$(OBJS): Makefile
-include ${DEPS}
