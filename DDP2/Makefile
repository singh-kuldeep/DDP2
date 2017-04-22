all: solver

solver: TVDMainSolver.cpp inputfile
	g++ TVDMainSolver.cpp -o solver -std=c++11 

.PHONY: clean
clean:
	rm -fr solver
