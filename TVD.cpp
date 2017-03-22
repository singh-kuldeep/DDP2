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
/// @image html Mach.eps
*/



/*! \file TVD.cpp
    \brief This header file contains the run() function which runs the solver. 
    \author Kuldeep Singh
    \date 2017
*/

#include "iostream"
#include <vector>
#include <fstream>
#include "math.h"
#include "time.h"
#include "netfluxinterface.h"
#include "BC.h"
#include "initial_condition.h"
// #include "local_time_step.h"

// Headers for grids 
#include "grid.h"
// #include "grid_ideal_nozzle.h"
// #include "grid_straight_duct.h"
// #include "grid_bump.h"
// #include "grid_nozzle.h"
// #include "grid_conical_nozzle.h"
// #include "grid_ideal_nozzle.h"

using namespace std ;

/** @brief This function implements the boundary condition, iFaceAreaVector is
	not required Because currently the flow in x direction and 2D flow */
void BC(
	vector<vector<vector<vector<double> > > > & ConservedVariables,
	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
	vector<vector<vector<vector<double> > > > & kFaceAreaVector,
	int Ni, int Nj, int Nk);

/** @brief This function generates the area vector and cell volumes inside the
	domain whole domain*/
void grid(
	vector<vector<vector<vector<double> > > > & iFaceAreaVector, 
	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
	vector<vector<vector<vector<double> > > > & kFaceAreaVector, 
	vector<vector<vector<double> > >& CellVolume,
	vector<vector<vector<double> > >& delta_s,
	int & Ni, int & Nj, int & Nk, int GeometryOption);

/*! \brief This function runs the solver.
    \warning Currently not using this, because  grid() is not calculating ds
    value properly. So recheck this function as well after fixing the grid()
    function.
	\return double
*/
int main()
{	
	/** \param GeometryOption Using this option the initial condition and the 
	grids(area vector and the cell volumes will be defined appropriately) */
	// int OptionGeometry = 1; // Straight duct 
	int OptionGeometry = 2; // Bump inside the straight duct
	// int OptionGeometry = 3; // Idel_Nozzle(Designed using MOC)
	// int OptionGeometry = 4; // Nozzle with basic initial condition

	time_t StartTime; /**\param StartTime Simulation starting time*/
	time_t EndTime ; /**\param EndTime Simulation ending time*/
	time(&StartTime); // noting the starting time

	double DeltaT = 0.000003; /**\param DeltaT Time step*/
	double TIME = 1e8*DeltaT;
	int IterationValues = 1e8; 
	/**\param IterationValues Total iterations = floor(TIME/DeltaT)*/

	int Ni;/**\param Ni Number of cells(Including ghosts) in in "i" direction.*/
	int Nj;/**\param Nj Number of cells(Including ghosts) in in "j" direction.*/
	int Nk;/**\param Nk Number of cells(Including ghosts) in in "k" direction.*/
	/** @brief Final value of the Ni,Nj,Nk has been decided inside the grid() 
	function. So, do not use these parameters untill the grid function is 
	callled*/

	// This will be while implementing the local time steeping
	#if 0
	double DeltaT = 0.0000015; // this is for CFL = 0.2
	int IterationValues = 1e4 ;
	
	double lenght = 26 ; // keep it even
	double delta = 1.0 ; // this basically defines the grid size

	int N = floor(lenght/delta);
	// extra 4 is added for ghost cell
	int Ni = (3*N + 4);
	int Nj = N+4;
	int Nk = 1+4; 
	// Because this is 2D-simulation so no need to take large number of grids 
	// in z direction 
	#endif  

	/**\param &iFaceAreaVector This is a pointer to the 4D vector which has the
	area vector of all faces which are in "i" direction.*/
	vector<vector<vector<vector<double> > > > iFaceAreaVector;
	/**\param &jFaceAreaVector This is a pointer to the 4D vector which has the
	area vector of all faces which are in "j" direction.*/
	vector<vector<vector<vector<double> > > > jFaceAreaVector;
	/**\param &kFaceAreaVector This is a pointer to the 4D vector which has the
	area vector of all faces which are in "k" direction.*/
	vector<vector<vector<vector<double> > > > kFaceAreaVector;
	/**\param CellVolumeIn Input pointer to cell volumes*/
	vector<vector<vector<double> > > CellVolume ;
	/**\param delta_s Minimum distance*/
	vector<vector<vector<double> > > delta_s ;
	
	grid(
		iFaceAreaVector,
		jFaceAreaVector,
		kFaceAreaVector,
		CellVolume,
		delta_s,
		Ni,Nj,Nk,OptionGeometry);

	// cout << "Ni, Nj, Nk :-> "<< Ni << "  " << Nj << "  " << Nk << endl;
	
	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	/**\param ConservedVariables This is the pointer to the 4D vector where all
	the conserved variables ([Density , x-momentum, y-momentum, z-momentum, 
	Energy]) of previous time step are stored.*/
	Dim4 ConservedVariables(Ni,Dim3(Nj,Dim2(Nk,Dim1(5)))); 

	/**\param ConservedVariablesNew This is the pointer to the 4D vector where
	all the conserved variables ([Density , x-momentum, y-momentum, z-momentum,
	Energy]) of current/new time step are stored.*/
	Dim4 ConservedVariablesNew(Ni,Dim3(Nj,Dim2(Nk,Dim1(5)))); 
	
	// Initializing the domain
	initial_condition(ConservedVariables, ConservedVariablesNew, Ni, Nj, Nk, OptionGeometry);
	
	#if 0 // for the testing purposes
	// storing the all conserved variables in one plane just after initialization
	ofstream kullu_2D_initial ;
	kullu_2D_initial.open("2D_parameters_B.csv");
	// kullu_2D_initial << "density" << "," << "density*u" << ","<< "density*v"
	// << "," << "density*w" << "," << "energy"  << endl ;
	for (int i = 2; i < Ni-2; ++i)
	{
		for (int j = 2; j < Nj-2; ++j)
		{
			kullu_2D_initial << ConservedVariables[i][j][Nk/2][0] << "," << 
			ConservedVariables[i][j][Nk/2][1] <<","<< 
			ConservedVariables[i][j][Nk/2][2] << "," <<
			ConservedVariables[i][j][Nk/2][3] << "," <<
			ConservedVariables[i][j][Nk/2][4] << endl ;
		}
	}
	return 0;
	#endif
	
	ofstream kullu_mass ;
	kullu_mass.open("Residual.csv");
	// kullu_mass <<  "t(secs)" << "," << "DensityResidual"  << "," <<
	// "xMomentumResidual" << "," << "yMomentumResidual" <<","<< 
	// "zMomentumResidual" << "," << "EnergyResidual" << endl ;

	double iCellInterfaceVolume;
	/**\param iCellInterfaceVolume Average of 
	right and left cell volume in i direction*/
	double jCellInterfaceVolume; 
	/**\param jCellInterfaceVolume Average of 
	right and left cell volume in j direction*/
	double kCellInterfaceVolume; 
	/**\param kCellInterfaceVolume Average of 
	right and left cell volume in k direction*/
	
	// Iterations starts here 
	// This file is opened to store the residuals at each time step
	for (int t = 0; t < IterationValues; ++t)
	{
		// cout << "timestep  = " << t << endl ; 
		
		// Before every time step we need to have proper value in the ghost 
		// cells So, BC takes  care of Inlet, Exit, y-wall and Z-wall BC
		BC(ConservedVariables,jFaceAreaVector,kFaceAreaVector,Ni,Nj,Nk); 
		
		// next time step ConservedVariabless calculation
		for (int i = 1; i < Ni-2; ++i)
		// Or (int i = 2; i < (Ni+1)-2; ++i) Total Ni+1 interface, 2 used in BC
		// implementation
		{
			for (int j = 1; j < Nj-2; ++j)
			{
				for (int k = 1; k < Nk-2; ++k)
				{	
					/**\bug Local time step needs to be used to reduce the 
					simulation time*/

					// DeltaT = dt(i+1,j+1,2,delta_s,ConservedVariables);
					// cout << "dt   " << DeltaT << endl;

					iCellInterfaceVolume = 0.5*(CellVolume[i][j][k] + 
						CellVolume[i+1][j][k]);
					jCellInterfaceVolume = 0.5*(CellVolume[i][j][k] + 
						CellVolume[i][j+1][k]);
					kCellInterfaceVolume = 0.5*( CellVolume[i][j][k] + 
						CellVolume[i][j][k+1]);

					// net flux using the class netfluxinterface
					netfluxinterface xrightface(
						ConservedVariables[i-1][j][k],
						ConservedVariables[i][j][k],
						ConservedVariables[i+1][j][k],
						ConservedVariables[i+2][j][k],
						iFaceAreaVector[i-1][j][k],
						iFaceAreaVector[i][j][k],
						iFaceAreaVector[i+1][j][k],
						CellVolume[i-1][j][k],
						CellVolume[i][j][k],
						CellVolume[i+1][j][k],
						CellVolume[i+2][j][k],
						DeltaT);

					netfluxinterface yrightface(
						ConservedVariables[i][j-1][k],
						ConservedVariables[i][j][k],
						ConservedVariables[i][j+1][k],
						ConservedVariables[i][j+2][k],
						jFaceAreaVector[i][j-1][k],
						jFaceAreaVector[i][j][k],
						jFaceAreaVector[i][j+1][k],
						CellVolume[i][j-1][k],
						CellVolume[i][j][k],
						CellVolume[i][j+1][k],
						CellVolume[i][j+2][k],
						DeltaT) ;

					netfluxinterface zrightface(
						ConservedVariables[i][j][k-1],
						ConservedVariables[i][j][k],
						ConservedVariables[i][j][k+1],
						ConservedVariables[i][j][k+2],
						kFaceAreaVector[i][j][k-1],
						kFaceAreaVector[i][j][k],
						kFaceAreaVector[i][j][k+1],
						CellVolume[i][j][k-1],
						CellVolume[i][j][k],
						CellVolume[i][j][k+1],
						CellVolume[i][j][k+2],
						DeltaT) ;

					// updating the ConservedVariablesNew using flux at the 
					// right interfaces
					for (int l = 0; l < 5; ++l)
					{
						ConservedVariablesNew[i][j][k][l] -=(DeltaT/
							iCellInterfaceVolume)*(xrightface.NetFlux[l]);
						ConservedVariablesNew[i+1][j][k][l] +=(DeltaT/
							iCellInterfaceVolume)*(xrightface.NetFlux[l]);

						ConservedVariablesNew[i][j][k][l] -=(DeltaT/
							jCellInterfaceVolume)*(yrightface.NetFlux[l]);
						ConservedVariablesNew[i][j+1][k][l] +=(DeltaT/
							jCellInterfaceVolume)*(yrightface.NetFlux[l]);

						ConservedVariablesNew[i][j][k][l] -=(DeltaT/
							kCellInterfaceVolume)*(zrightface.NetFlux[l]);
						ConservedVariablesNew[i][j][k+1][l] +=(DeltaT
							/kCellInterfaceVolume)*(zrightface.NetFlux[l]);
					}
				}
			}
			
		}

		// Residual calculation after each time step and writing the all 
		//residuals into the file
		double DensityResidual = 0.0; 
		/**\param DensityResidual Density residual*/
		double xMomentumResidual = 0.0 ; 
		/**\param xMomentumResidual x Momentum residual*/
		double yMomentumResidual = 0.0 ; 
		/**\param yMomentumResidual y Momentum residual*/
		double zMomentumResidual = 0.0 ; 
		/**\param zMomentumResidual z Momentum residual*/
		double EnergyResidual = 0.0 ; 
		/**\param Energy residual */

		int TotalGridPoints = 0 ; 
		for (int x = 2; x < Ni-2; ++x)
		{
			for (int y = 2; y < Nj-2; ++y)
			{
				TotalGridPoints   += 1 ; 
				DensityResidual   += pow((ConservedVariablesNew[x][y][2][0] -
					ConservedVariables[x][y][2][0]),2);
				xMomentumResidual += pow((ConservedVariablesNew[x][y][2][1] - 
					ConservedVariables[x][y][2][1]),2);     
				yMomentumResidual += pow((ConservedVariablesNew[x][y][2][2] - 
					ConservedVariables[x][y][2][2]),2);     
				zMomentumResidual += pow((ConservedVariablesNew[x][y][2][3] - 
					ConservedVariables[x][y][2][3]),2);     
				EnergyResidual    += pow((ConservedVariablesNew[x][y][2][4] - 
					ConservedVariables[x][y][2][4]),2);     
			}
		}
		if (t%10 == 0)
		{
			// cout << "TotalGridPoints" << TotalGridPoints << endl ;
			kullu_mass << t << "," << t*DeltaT << "," << sqrt(DensityResidual/
			((Ni-4)*(Nj-4)))  << "," << sqrt(xMomentumResidual/((Ni-4)*(Nj-4)))
			<< "," << sqrt(yMomentumResidual/((Ni-4)*(Nj-4))) <<","<< sqrt(
			zMomentumResidual/((Ni-4)*(Nj-4))) << "," << sqrt(EnergyResidual/
			((Ni-4)*(Nj-4))) << endl ;			
		}

		if (t%10==0)
		{
			cout <<  t << "  --->  " << "  "<<  sqrt(DensityResidual) << endl ;
		}

		// before going to the new time step updating the old conserved 
		// variables by new ones.
		for (int i = 2; i < Ni-2; ++i)
		{
			for (int j = 2; j < Nj-2; ++j)
			{
				for (int k = 2; k < Nk-2; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						ConservedVariables[i][j][k][l] = 
						ConservedVariablesNew[i][j][k][l] ;
					}				
				}
			}
		}

		if (t%1 == 0)
		{
			// storing the all conserved variables in one plane
			ofstream kullu_2D ;
			kullu_2D.open("2D_parameters_B.csv");
			// kullu_2D << "density" << "," << "density*u" << ","<< "density*v"
			// << "," << "density*w" << "," << "energy"  << endl ;
			for (int i = 2; i < Ni-2; ++i)
			{
				for (int j = 2; j < Nj-2; ++j)
				{
					kullu_2D << ConservedVariables[i][j][Nk/2][0] << "," << 
					ConservedVariables[i][j][Nk/2][1] <<","<< 
					ConservedVariables[i][j][Nk/2][2] << "," <<
					ConservedVariables[i][j][Nk/2][3] << "," <<
					ConservedVariables[i][j][Nk/2][4] << endl ;
				}
			}
		}

		// Restart file is being written to restart the simulation where it was
		// left earlier
		if (t%500 == 0)
		{
			// storing the conserved variables throughout the domain
			ofstream kullu_restart ;
			kullu_restart.open("restart.csv");
			// kullu_restart << "density" << "," << "density*u" << ","<<
			// "density*v" << "," << "density*w" << "," << "energy"  << endl ;
			for (int i = 0; i < Ni; ++i)
			{
				for (int j  = 0; j < Nj; ++j)
				{
					for (int k = 0; k < Nk; ++k)
					{
						kullu_restart << ConservedVariables[i][j][k][0] << endl
						<< ConservedVariables[i][j][k][1] << endl<< 
						ConservedVariables[i][j][k][2] << endl << 
						ConservedVariables[i][j][k][3] << endl << 
						ConservedVariables[i][j][k][4] << endl ;
					}
				}
			}
		}
	} 
	// time progression ends here 
	time(&EndTime) ;
	double SimulationTotalTime = difftime (EndTime,StartTime);
	cout << "Time taken by the solver in secs = " << SimulationTotalTime<<endl;
	// return 0;
	return 0;
}
