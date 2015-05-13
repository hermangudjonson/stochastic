
# Makefile for Branching Process Simulation Project
# homebrew c++ compiler (in /usr/local/Cellar/gcc/4.9.1)
CC=g++-4.9

# compiler flags
CFLAGS=-W -Wno-long-long -pedantic -Wno-variadic-macros -std=c++11

INCLUDE=

all: basic_gillespie.o target

basic_gillespie.o: basic_gillespie.cpp basic_gillespie.h timekeeper.h
	$(CC) $(CFLAGS) $(INCLUDE) -c basic_gillespie.cpp

target: yan_network.cpp basic_gillespie.o
	$(CC) $(CFLAGS) $(INCLUDE) yan_network.cpp basic_gillespie.o -o yan
clean: 
	rm -rf *.o testing *~
