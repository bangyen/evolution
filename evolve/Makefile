# Makefile to compile seria evolution code
#       Joshua Hoskins 
#

LIB        = -L/usr/lib64/ -L/usr/lib/
INCLUDES   = -Iinclude/ -I/usr/include/
CC         = g++
SRC        = src
CFLAGS     = -O -std=c++11 -Wall

all: evolution

%.o: %.cc
	${CC} ${CFLAGS} -c -o $@ $< 
evolution : plot.cpp ${SRC}/evolve-serial.o ${SRC}/evolve.o ${SRC}/runge-kutta.o ${SRC}/gpd.cpp
	${CC} ${INCLUDES} -o $@  ${CFLAGS} $^ ${LIB}
clean:
	rm -f *.o *~
