/*! \mainpage My Personal Index Page
*
*\section intro_sec Introduction
*
* This is the C++ code to solve the high speed fluid flow. Currently, Euler 
* flow is being solved but this code has been designed in moulder way so to  
* solve the viscus flow additional viscus flux class can be added very easily. 
* This code has been written to fulfill the requirement of the Dual Degree
* Project(DDP). 
*
*\section install_sec Installation & Use
*
* To use the solver. Follow these simple steps.
*  - Download form here : https://github.com/singh-kuldeep/DDP2 or 
* 	<a href="https://github.com/singh-kuldeep/DDP2">click here</a>
*  - Go to the folder DDP2 and compile and run the file TVD.cpp 
* 	(ex. g++ TVD.cpp && ./a.out)  
*  - Nozzle has been set up as a default geometry but it can be changed from 
*	"run.h" file by uncommenting the header file  
*  - Currently there are two different geometry options are available
*     1. Curved wall high area ratio diverging nozzle
*     2. Triangular bump inside straight duct 
* 
*\section brief Brief About the solver 
*  - 3D 
*  - Roe scheme based
*  - C++
*  - Exact theory can be found <a href="https://drive.google.com/open?id=0B9x_nh0D_HhzMnBjc0w5MmJpcnc">here</a>  
*
*\section input Input to the solver 
*  - Grid points
*  - Boundary condition
*  - Some initial condition
*/

#include <iostream>
#include "run.h"
using namespace std;
int main()
{
	run();
	return 0;
}
