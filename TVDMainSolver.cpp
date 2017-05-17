/*! \mainpage My Personal Index Page

\section intro_sec Introduction

 This is the C++ code to solve the high speed fluid flow. Currently, Euler 
 flow is being solved but this code has been designed in moulder way so to  
 solve the viscus flow additional "viscus flux" class can be added very easily. 
 This code has been written to fulfill the requirement of the Dual Degree
 Project(DDP). All the major input has been given from the input file. For post 
 processing either MATLAB or Python could be used.  

\section install_sec Installation & Use

 To use the solver. Follow these simple steps.
	- Download form here : https://github.com/singh-kuldeep/DDP2 or 
 	<a href="https://github.com/singh-kuldeep/DDP2">click here</a>
	- Go to the folder DDP2 and compile and run the file TVD.cpp 
 	(ex. g++ MainSolver.cpp -std=c++11 && ./a.out)  
	- Nozzle has been set up as a default geometry but it can be changed from 
	"inputfile" by uncommenting the appropriate geometry  

\section input Input file
All the major inputs are taken through the input file. 
	1. CFL (Courant–Friedrichs–Lewy)
	2. Currently there are three different geometry options are available
     - CD Nozzle generated using MOC(Default)
     - Triangular bump inside straight duct
     - Straight rectangular duct
 
\section brief Brief about the solver 
	- 3D Cartesian (x,y,z)  
	- Roe scheme based
	- C++
	- Exact theory can be found 
<a href="https://drive.google.com/open?id=0B9x_nh0D_HhzMnBjc0w5MmJpcnc">here</a>  

\section input Input to the solver 
	- Grid points
	- Boundary condition
	- Some initial condition
\section output Output files.
Here are the list of files which will come as the output of the solver.
	- Residual_Nozzle.csv	: This file contains the all the residuals 
(Mass, Momentum, Energy).
	- grids_Nozzle_2D.csv	: This file contains the grid point 
(x,y) coordinates. 
	- 2D_parameters_B.csv	: This file contains all the conserved parameters 
at the 2D plane. 
\section plot Results & Plots
Same older contains the MATLAB script "plot_data.m". Once the simulation has 
started and the output files are 
generated, one can simply run the MATALB script and can see the plots which 
are listed below.   
	- Density Residual
	- X Momentum Residual 
	- Y Momentum Residual 
	- Z Momentum Residual
	- Energy Residual 
	- Mach Number 
	- Density
	- Velocity
	- Temperature
	- Pressure 
	- Geometry 2D cross section	
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

#include "flux.h"
#include "boundaryNetflux.h"
#include "deltat.h"
#include "residual.h"
#include "initial_condition.h"
#include "grid.h" // Headers for grids 
#include "ghostcell.h" // Headers for ghost cells
#include "array_tester.h" // test the 3D/4D array
#include "WriteConservedQuantities.h"
#include "colortext.h" 
using namespace std ;


int main()
{		
	time_t StartTime; /**\param StartTime Simulation starting time*/
	time_t EndTime ; /**\param EndTime Simulation ending time*/
	time(&StartTime); // noting the starting time

	/**\param TotalIteration Total iterations */
	int TotalIteration;
	
	string Scheme ;// = "AUSM";// or "Roe"/
	string TimeSteping ;
	string GeometryOption ;
	string gamma; 
	/** \param GeometryOption Using this option grids (area vector and the 
	cell volumes) will be defined appropriately */
	
	double CFL; /**\param CFL */
	double SpecificHeatRatio;

	/**\param deltat Time step */
	double deltat = 1000; 

	ifstream infile("inputfile");
	string aline;
	// reading the input file
	while(!infile.eof())// file ended
	{
		getline(infile,aline); // reading line form file

		if (aline.find( "//" )!=0 && aline.empty()==false) 
		{
			if(aline.find("TotalIteration")!=string::npos)
			{
				TotalIteration = stoi (aline.substr(aline.find("=")+1));
			}
			else if (aline.find("Scheme")!=string::npos)
			{
				Scheme = aline.substr(aline.find("=")+2); 
			}
			else if (aline.find("GeometryOption")!=string::npos)
			{
				GeometryOption = aline.substr(aline.find("=")+2);
			}				
			else if(aline.find("CFL")!=string::npos)
			{
				CFL = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("CFL")!=string::npos)
			{
				CFL = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("TimeSteping")!=string::npos)
			{
				TimeSteping = aline.substr(aline.find("=")+2);
			}
			else if(aline.find("gamma")!=string::npos)
			{
				gamma = aline.substr(aline.find("=")+2);
			}
			else if (aline.find("SpecificHeatRatio")!=string::npos)
			{
				SpecificHeatRatio = stod (aline.substr(aline.find("=")+1)); 
			}								
		}
	}
	// Reading the input file over


	int Ni;/**\param Ni Number of live cells in in "i" direction.*/
	int Nj;/**\param Nj Number of live cells in in "j" direction.*/
	int Nk;/**\param Nk Number of live cells in in "k" direction.*/
	/** @brief Final value of the Ni,Nj,Nk has been decided inside the grid() 
	function. So, do not use these parameters until the grid function is 
	called*/

	/**\param Coordinates This is a 4D vector which has the
	coordinated of all vertices */
	vector<vector<vector<vector<double> > > > Coordinates;
	/**\param iFaceAreaVector This is a 4D vector which has the
	area vector of all faces which are in "i" direction.*/
	vector<vector<vector<vector<double> > > > iFaceAreaVector;
	/**\param jFaceAreaVector This is a 4D vector which has the
	area vector of all faces which are in "j" direction.*/
	vector<vector<vector<vector<double> > > > jFaceAreaVector;
	/**\param kFaceAreaVector This is a 4D vector which has the
	area vector of all faces which are in "k" direction.*/
	vector<vector<vector<vector<double> > > > kFaceAreaVector;
	/**\param CellVolume Input pointer to cell volumes*/
	vector<vector<vector<double> > > CellVolume ;
	/**\param delta_s Minimum distance*/
	vector<vector<vector<double> > > delta_s ;
	// After calling grid function all the live cell quantities will be decided
	grid(
		Coordinates,
		iFaceAreaVector,
		jFaceAreaVector,
		kFaceAreaVector,
		CellVolume,
		delta_s,
		Ni,Nj,Nk,
		GeometryOption);

	// if(test4Darray("iFaceAreaVector",iFaceAreaVector,Ni+1,Nj,Nk,3)==0)
	// {
	// 	return 0;
	// }
	// After this point the Ni, Nj, Nk has been decided 
	cout << "Ni, Nj, Nk :-> "<< Ni << "  " << Nj << "  " << Nk << endl;
	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	// checking the NaN in the grid parameters
	#if 0
	if(test4DArray("Coordinates",Coordinates,Ni,Nj,Nk,3) == 0)
	{
		return 0;
	}
 	if(test4DArray("iFaceAreaVector",iFaceAreaVector,Ni+1,Nj,Nk,3)==0)
 	{
 		return 0;
 	} 
 	if(test4DArray("jFaceAreaVector",jFaceAreaVector,Ni,Nj+1,Nk,3)==0)
 	{
 		return 0;
 	} 
 	if(test4DArray("kFaceAreaVector",kFaceAreaVector,Ni,Nj,Nk+1,3)==0)
 	{
 		return 0;
 	} 

 	if(test3DArray("CellVolume",CellVolume,Ni,Nj,Nk)==0)
 	{
 		return 0;
 	}	
 	if(test3DArray("delta_s",delta_s,Ni,Nj,Nk)==0)
 	{
 		return 0;
 	}
 	#endif
	/**\param i0GhostCellVolume Ghost cell volume array at i = 0*/
	/**\param j0GhostCellVolume Ghost cell volume array at j = 0*/
	/**\param k0GhostCellVolume Ghost cell volume array at k = 0*/
	/**\param iNiGhostCellVolume Ghost cell volume array at i = Ni*/
	/**\param jNjGhostCellVolume Ghost cell volume array at j = Nj*/
	/**\param kNkGhostCellVolume Ghost cell volume array at k = Nk*/
	
	Dim3 i0GhostCellVolume(1,Dim2(Nj,Dim1(Nk))); 
	Dim3 iNiGhostCellVolume(1,Dim2(Nj,Dim1(Nk))); 
	Dim3 j0GhostCellVolume(Ni,Dim2(1,Dim1(Nk))); 
	Dim3 jNjGhostCellVolume(Ni,Dim2(1,Dim1(Nk))); 
	Dim3 k0GhostCellVolume(Ni,Dim2(Nj,Dim1(1))); 
	Dim3 kNkGhostCellVolume(Ni,Dim2(Nj,Dim1(1))); 
	
	/*After calling ghostcell function volumes of all the ghost cells 
	will be decided*/
	ghostcell(Coordinates,iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,
	CellVolume, i0GhostCellVolume,j0GhostCellVolume, k0GhostCellVolume,
	iNiGhostCellVolume,jNjGhostCellVolume,kNkGhostCellVolume, Ni, Nj, Nk);

	#if 0
	// checking the NaN 
 	if(test3DArray("i0GhostCellVolume",i0GhostCellVolume,1,Nj,Nk)==0)
 	{
 		return 0;
 	}	
 	if(test3DArray("iNiGhostCellVolume",iNiGhostCellVolume,1,Nj,Nk)==0)
 	{
 		return 0;
 	}
 	if(test3DArray("j0GhostCellVolume",j0GhostCellVolume,Ni,1,Nk)==0)
 	{
 		return 0;
 	}
 	if(test3DArray("jNjGhostCellVolume",jNjGhostCellVolume,Ni,1,Nk)==0)
 	{
 		return 0;
 	}
 	if(test3DArray("k0GhostCellVolume",k0GhostCellVolume,Ni,Nj,1)==0)
 	{
 		return 0;
 	}
 	if(test3DArray("kNkGhostCellVolume",kNkGhostCellVolume,Ni,Nj,1)==0)
 	{
 		return 0;
 	}
	#endif
	/**\param ConservedVariables This is the pointer to the 4D vector where all
	the conserved variables ([Density , x-momentum, y-momentum, z-momentum, 
	Energy]) of previous time step are stored.*/
	Dim4 ConservedVariables(Ni,Dim3(Nj,Dim2(Nk,Dim1(5)))); 

	/**\param ConservedVariablesNew This is the pointer to the 4D vector where
	all the conserved variables ([Density , x-momentum, y-momentum, z-momentum,
	Energy]) of current/new time step are stored.*/
	Dim4 ConservedVariablesNew(Ni,Dim3(Nj,Dim2(Nk,Dim1(5)))); 
	
	// Initializing the domain
	initial_condition(ConservedVariables, ConservedVariablesNew,
	 Ni, Nj, Nk, SpecificHeatRatio);
	

	#if 0
	// checking whether initializations is proper  
	if(testConservedVariables("ConservedVariables",ConservedVariables,Ni,Nj,Nk,5) == 0)
	{
		return 0;
	}
	if(testConservedVariables("ConservedVariablesNew",ConservedVariablesNew,Ni,Nj,Nk,5) == 0)
	{
		return 0;
	}
	#endif

	/**\param iFacesFlux This is a 4D vector which has the
	fluxes of all faces which are in "i" direction.*/
	Dim4 iFacesFlux(Ni+1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	/**\param jFacesFlux This is a 4D vector which has the
	fluxes of all faces which are in "j" direction.*/
	Dim4 jFacesFlux(Ni,Dim3(Nj+1,Dim2(Nk,Dim1(5))));
	/**\param kFacesFlux This is a 4D vector which has the
	fluxes of all faces which are in "k" direction.*/
	Dim4 kFacesFlux(Ni,Dim3(Nj,Dim2(Nk+1,Dim1(5))));

	ofstream kullu_Residual ;
	kullu_Residual.open("./Results/outputfiles/Residual.csv");
	// kullu_Residual <<  "t(secs)" << "," << "DensityResidual"  << "," <<
	// "xMomentumResidual" << "," << "yMomentumResidual" <<","<< 
	// "zMomentumResidual" << "," << "EnergyResidual" << endl ;

	//structure is i=0,i=Ni,j=0,j=Nj,k=0,k=Nk
	ofstream Mass_Residual ;
	Mass_Residual.open("./Results/outputfiles/MassResidualAllBoundary.csv");
	ofstream XMom_Residual ;
	XMom_Residual.open("./Results/outputfiles/XMomResidualAllBoundary.csv");
	ofstream YMom_Residual ;
	YMom_Residual.open("./Results/outputfiles/YMomResidualAllBoundary.csv");
	ofstream Zmom_Residual ;
	Zmom_Residual.open("./Results/outputfiles/ZmomResidualAllBoundary.csv");
	ofstream Energy_Residual ;
	Energy_Residual.open("./Results/outputfiles/EnergyResidualAllBoundary.csv");
	
	// Iterations starts here 
	for (int iteration = 0; iteration < TotalIteration; ++iteration)
	{	
		double deltat;
		// calculating the global time step after every time iteration
		if(TimeSteping == "Global")
		{
			deltat = getGlobalDeltaT(ConservedVariables, delta_s, CFL, Ni, Nj, 
				Nk, gamma, SpecificHeatRatio);
		}
		
		// flux
		flux(iFacesFlux,jFacesFlux,kFacesFlux,iFaceAreaVector,jFaceAreaVector,
			kFaceAreaVector,ConservedVariables,Ni,Nj,Nk,gamma,SpecificHeatRatio);

		// return 0;





		// updating the conserved variables 
		for (int i = 0; i < Ni; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						// Local time step
						if(TimeSteping == "Local")
						{
							deltat = getLocalDeltaT(ConservedVariables[i][j][k],
								delta_s[i][j][k],CFL,gamma,SpecificHeatRatio);
						}
						vector<double> NetFlux(5);
						NetFlux[l] = (iFacesFlux[i][j][k][l] - iFacesFlux[i+1][j][k][l]) + 
						(jFacesFlux[i][j][k][l] - jFacesFlux[i][j+1][k][l]) + 
						(kFacesFlux[i][j][k][l] - kFacesFlux[i][j][k+1][l]);
						
						ConservedVariablesNew[i][j][k][l] +=(deltat/CellVolume[i][j][k])*
						NetFlux[l]; 
					}
				}
			}
		}

		#if 0
		if(testConservedVariables("ConservedVariables",ConservedVariables,Ni,Nj,Nk,5) == 0)
		{
			cout << "Check the line number "<< __LINE__ << " in TVDMainSolver.cpp" << endl;
			return 0;
		}
		if(testConservedVariables("ConservedVariablesNew",ConservedVariablesNew,Ni,Nj,Nk,5) == 0)
		{
			cout << "Check the line number "<< __LINE__ << " in TVDMainSolver.cpp" << endl;
			return 0;
		}
		

		#endif
		// if(iteration%10 == 0)
		// {
		// 	cout << "iteration -> " << iteration << " dt -> " << deltat; 
		// }
			
		if(iteration%10 == 0)
		{
			cout << iteration << "  " << deltat << "  " ; 
			std::vector<double> Residual(5);
			residual(Residual,iteration,ConservedVariables,ConservedVariablesNew,Ni,Nj,Nk);
			kullu_Residual << iteration << "," << Residual[0] << "," << Residual[1]<<"," << 
			Residual[2] <<"," << Residual[3] <<","<< Residual[4]<< ","<< Residual[5]<< endl;
		}
		
		#if 1
		if(iteration%10 == 0)
		{
			// Net Fluxes integration at theboundaries 
			vector<double> iNetFlux(5);
			vector<double> iNiNetFlux(5);
			vector<double> jNetFlux(5);
			vector<double> jNjNetFlux(5);
			vector<double> kNetFlux(5);
			vector<double> kNkNetFlux(5);

			boundaryNetflux(iNetFlux, iNiNetFlux, jNetFlux, jNjNetFlux, 
				kNetFlux, kNkNetFlux, iFacesFlux,jFacesFlux,kFacesFlux,
			iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,Ni, Nj, Nk);

			Mass_Residual << iNetFlux[0]<<","<<iNiNetFlux[0]<<","<<jNetFlux[0]<<
			","<<jNjNetFlux[0]<<","<<kNetFlux[0]<<","<<kNkNetFlux[0] << endl;
			XMom_Residual << iNetFlux[1]<<","<<iNiNetFlux[1]<<","<<jNetFlux[1]<<
			","<<jNjNetFlux[1]<<","<<kNetFlux[1]<<","<<kNkNetFlux[1] << endl;
			YMom_Residual << iNetFlux[2]<<","<<iNiNetFlux[2]<<","<<jNetFlux[2]<<
			","<<jNjNetFlux[2]<<","<<kNetFlux[2]<<","<<kNkNetFlux[2] << endl;
			Zmom_Residual << iNetFlux[3]<<","<<iNiNetFlux[3]<<","<<jNetFlux[3]<<
			","<<jNjNetFlux[3]<<","<<kNetFlux[3]<<","<<kNkNetFlux[3] << endl;
			Energy_Residual << iNetFlux[4]<<","<<iNiNetFlux[4]<<","<<jNetFlux[4]<<
			","<<jNjNetFlux[4]<<","<<kNetFlux[4]<<","<<kNkNetFlux[4] << endl;	
		}
		#endif

		// before going to the new time step updating the old conserved 
		// variables by new ones.
		for (int i = 0; i < Ni; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						ConservedVariables[i][j][k][l] = 
						ConservedVariablesNew[i][j][k][l] ;
					}				
				}
			}
		}

		WriteConserveredQuantities(ConservedVariables,Ni,Nj,Nk);
	} 

	// time progression ends here 
	time(&EndTime) ;
	double SimulationTotalTime = difftime (EndTime,StartTime);
	cout << "Time taken by the solver in secs = " << SimulationTotalTime<<endl;
	// return 0;
	return 0;
}
