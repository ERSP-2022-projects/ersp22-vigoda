default: main

main: treegen

treegen: treegen.cpp tree.o NewickTree.o
	g++ --std=c++17 treegen.cpp tree.o NewickTree.o -o treegen.out
	mkdir -p results

NewickTree.o: NewickTree.cpp NewickTree.h
	g++ --std=c++17 -c NewickTree.cpp 

generateMixture: generateMixture.cpp NewickTree.o
	g++ --std=c++17 generateMixture.cpp NewickTree.o -o generateMixture.out

newicktest: newicktest.cpp NewickTree.o
	g++ --std=c++17 newicktest.cpp NewickTree.o -o newicktest

tree.o: tree.h tree.cpp
	g++ --std=c++17 -c tree.cpp

clean: 
	rm -f *.o
	rm -f *.out
	rm -r results/*
