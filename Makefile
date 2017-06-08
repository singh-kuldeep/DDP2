all: solver

solver: MainSolver.cpp inputfile
	g++ MainSolver.cpp -o solver  

.PHONY: clean
clean:
	rm -fr solver
