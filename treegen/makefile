default: treegen

all: treegen generateMixture

treegen: treegen.cpp tree.o NewickTree.o
	g++ --std=c++17 -I "C:\Eigen" -g treegen.cpp tree.o NewickTree.o -o treegen.out
	mkdir -p results

NewickTree.o: NewickTree.cpp NewickTree.h
	g++ --std=c++17 -c -I "C:\Eigen" NewickTree.cpp 

generateMixture: generateMixture.cpp NewickTree.o
	g++ --std=c++17 -I "C:\Eigen" generateMixture.cpp NewickTree.o -o generateMixture.out

newicktest: newicktest.cpp NewickTree.o
	g++ --std=c++17 -g newicktest.cpp NewickTree.o -o newicktest

tree.o: tree.h tree.cpp
	g++ --std=c++17 -c -g tree.cpp

EigenTest: EigenTest.cpp
	g++ --std=c++11 -I "C:\Eigen" EigenTest.cpp -o EigenTest.out

debugGenerateMixture: generateMixture.cpp NewickTree.o
	g++ --std=c++17 -g generateMixture.cpp NewickTree.o -o debugGenerateMixture.out

varyBranchLength: varyBranchLength.cpp NewickTree.o
	g++ --std=c++17 -O3 -I "C:\Eigen" varyBranchLength.cpp NewickTree.o -o varyBranchLength.out

clean: 
	rm -f *.o
	rm -f *.out

clean_results:
	rm -r results/*
clean_all: clean clean_results
debug_treegen: 
	gdb treegen.out
debug_mixture: 
	gdb generateMixture.out
