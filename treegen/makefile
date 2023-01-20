default: main

main: main.cpp NewickTree.o
	g++ --std=c++11 main.cpp NewickTree.o -o main.out
	mkdir -p results

newicktest: newicktest.o NewickTree.o
	g++ --std=c++11 newicktest.o NewickTree.o -o newicktest
newicktest.o: newicktest.cpp NewickTree.h
	g++ --std=c++11 -c newicktest.cpp
NewickTree.o: NewickTree.cpp NewickTree.h
	g++ --std=c++11 -c NewickTree.cpp
clean: 
	rm -f *.o
	rm -f *.out
	rm -f results/*.txt
	rm -f results/*.nex