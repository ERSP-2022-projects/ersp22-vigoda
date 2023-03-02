default: main

main: treegen

treegen: treegen.cpp tree.o NewickTree.o
	g++ --std=c++17 treegen.cpp tree.o NewickTree.o -o treegen.out
	mkdir -p results

NewickTree.o: NewickTree.cpp NewickTree.h
	g++ --std=c++17 -c -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 NewickTree.cpp 

generateMixture: generateMixture.cpp NewickTree.o
	g++ --std=c++17 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 generateMixture.cpp NewickTree.o -o generateMixture.out

newicktest: newicktest.cpp NewickTree.o
	g++ --std=c++17 newicktest.cpp NewickTree.o -o newicktest

tree.o: tree.h tree.cpp
	g++ --std=c++17 -c tree.cpp

EigenTest: EigenTest.cpp
	g++ --std=c++11 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 EigenTest.cpp -o EigenTest.out

debugGenerateMixture: generateMixture.cpp NewickTree.o
	g++ --std=c++17 -g generateMixture.cpp NewickTree.o -o debugGenerateMixture.out

varyBranchLength: varyBranchLength.cpp NewickTree.o
	g++ --std=c++17 -O3 -I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 varyBranchLength.cpp NewickTree.o -o varyBranchLength.out

clean: 
	rm -f *.o
	rm -f *.out
