main: main.cpp generator.o
	g++ main.cpp generator.o -o main.out

generator.o: generator.h generator.cpp rng.h consts.h
	g++ generator.cpp -c

clean: 
	rm -f *.o
	rm -f *.out
#	rm -f *.txt