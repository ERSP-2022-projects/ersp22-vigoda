default: all

all: treegen mixture 

main: treegen

treegen: treegen.cpp tree.o
	g++ --std=c++17 -g -Wall treegen.cpp tree.o -o treegen.out
	mkdir -p results

mixture: mixture.cpp tree.o
	g++ --std=c++17 -g -Wall mixture.cpp tree.o -o mixture.out

tree.o: tree.h tree.cpp
	g++ --std=c++17 -g -Wall -c tree.cpp

clean: 
	rm -f *.o
	rm -f *.out
	rm -r results/*
