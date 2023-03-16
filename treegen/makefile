default: clean

main: treegen

treegen: treegen.cpp tree.o
	g++ --std=c++17 treegen.cpp tree.o -o treegen.out
	mkdir -p results

tree.o: tree.h tree.cpp
	g++ --std=c++17 -c tree.cpp

resultreader: resultreader.cpp
	g++ --std=c++17 resultreader.cpp -o resultreader.out
	mkdir -p sumt

clean: 
	rm -f *.o
	rm -f *.out
	rm -r results/*
