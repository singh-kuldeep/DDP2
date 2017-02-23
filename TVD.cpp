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
*	- Download form here : https://github.com/singh-kuldeep/DDP2 or 
* 	<a href="https://github.com/singh-kuldeep/DDP2">click here</a>
*	- Go to the folder DDP2 and compile and run the file TVD.cpp 
* 	(ex. g++ TVD.cpp && ./a.out)  
*	- Nozzle has been set up as a default geometry but it can be changed from 
*	"run.h" file by uncommenting the header file  
*	- Currently there are two different geometry options are available
*     1. Curved wall high area ratio diverging nozzle
*     2. Triangular bump inside straight duct 
* 
*\section brief Brief about the solver 
*	- 3D Cartesian (x,y,z)  
*	- Roe scheme based
*	- C++
*	- Exact theory can be found 
<a href="https://drive.google.com/open?id=0B9x_nh0D_HhzMnBjc0w5MmJpcnc">here</a>  
*
*\section input Input to the solver 
*	- Grid points
*	- Boundary condition
*	- Some initial condition
*\section output Output files.
*Here are the list of files which will come as the output of the solver.
*	- Residual_Nozzle.csv	: This file contains the all the residuals 
*(Mass, Momentum, Energy).
*	- grids_Nozzle_2D.csv	: This file contains the grid point 
*(x,y) coordinates. 
*	- 2D_parameters_B.csv	: This file contains all the conserved parameters 
*at the 2D plane. 
*\section plot Results & Plots
*Same older contains the MATLAB script "plot_data.m". Once the simulation has 
*started and the output files are 
*generated, one can simply run the MATALB script and can see the plots which 
*are listed below.   
*	- Density Residual
*	- X Momentum Residual 
*	- Y Momentum Residual 
*	- Z Momentum Residual
*	- Energy Residual 
*	- Mach Number 
*	- Density
*	- Velocity
*	- Temperature
*	- Pressure 
*	- Geometry 2D cross section	
*/

#include "run.h"
int main()
{
	run();
	return 0;
}
