all: program

#program:	integral.o
#	g++ -Wall integral.o -o program.exe
	
program: integral.cpp	
	g++ -Wall -std=c++17 -fopenmp -o integral ./integral.cpp
		
clean:
	rm -f *.o program.exe