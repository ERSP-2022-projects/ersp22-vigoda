default: treegen

treegen: treegen.cpp NewickTree.o
	g++ --std=c++17 treegen.cpp NewickTree.o -o treegen.out
	mkdir -p results

NewickTree.o: NewickTree.cpp NewickTree.h
	g++ --std=c++17 -c NewickTree.cpp 
clean: 
	rm -f *.o
	rm -f *.out
	rm -f results/*.txt
	rm -f results/*.nex