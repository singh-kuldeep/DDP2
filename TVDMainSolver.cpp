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
#include "netfluxAUSM.h" // AUSM
#include "netfluxRoe.h" // Roe
#include "BC.h"
#include "initial_condition.h"

#include "grid.h" // Headers for grids 
#include "ghostcell.h" // Headers for ghost cells

using namespace std ;

// /** @brief This function implements the boundary condition, iFaceAreaVector is
// 	not required Because currently the flow in x direction and 2D flow */
// void BC(
// 	vector<vector<vector<vector<double> > > > & ConservedVariables,
// 	vector<vector<vector<vector<double> > > > & iFaceAreaVector,
// 	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
// 	vector<vector<vector<vector<double> > > > & kFaceAreaVector,
// 	int Ni, int Nj, int Nk, string InletBC, double TemperatureFreestream,
// 	double PressureFreestream, double MachFreestream);

/** @brief This function generates the area vector and cell volumes inside the
	domain whole domain*/
// void grid(
// 	vector<vector<vector<vector<double> > > > & Coordinates, 
// 	vector<vector<vector<vector<double> > > > & iFaceAreaVector, 
// 	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
// 	vector<vector<vector<vector<double> > > > & kFaceAreaVector, 
// 	vector<vector<vector<double> > >& CellVolume,
// 	vector<vector<vector<double> > >& delta_s,
// 	int & Ni, int & Nj, int & Nk, string GeometryOption);

// void ghostcell(
// 	vector<vector<vector<vector<double> > > > Coordinates,
// 	vector<vector<vector<vector<double> > > > iFaceAreaVector,
// 	vector<vector<vector<vector<double> > > > jFaceAreaVector,
// 	vector<vector<vector<vector<double> > > > kFaceAreaVector,
// 	vector<vector<vector<double> > > CellVolume,

// 	vector<vector<vector<vector<double> > > > & i0GhostCellVolume,
// 	vector<vector<vector<vector<double> > > > & j0GhostCellVolume,
// 	vector<vector<vector<vector<double> > > > & k0GhostCellVolume,

// 	vector<vector<vector<vector<double> > > > & iNiGhostCellVolume,
// 	vector<vector<vector<vector<double> > > > & jNjGhostCellVolume,
// 	vector<vector<vector<vector<double> > > > & kNkGhostCellVolume,
// 	int Ni, int Nj, int Nk)

/*! \brief This function runs the solver.
    \warning Currently not using this, because  grid() is not calculating ds
    value properly. So recheck this function as well after fixing the grid()
    function.
	\return double
*/
int main()
{		
	time_t StartTime; /**\param StartTime Simulation starting time*/
	time_t EndTime ; /**\param EndTime Simulation ending time*/
	time(&StartTime); // noting the starting time

	/**\param TotalIteration Total iterations = floor(TIME/deltat)*/
	int TotalIteration;
	
	string Scheme ;// = "AUSM";// or "Roe"/
	string InletBC ;
	string ExitBC ;
	string InitialCondition ;
	string TimeSteping ;
	string GeometryOption ;
	/** \param GeometryOption Using this option grids (area vector and the 
	cell volumes) will be defined appropriately */
	
	double CFL; /**\param CFL */

	// double TIME = TotalIteration*deltat;
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
			else if(aline.find("TimeSteping")!=string::npos)
			{
				TimeSteping = aline.substr(aline.find("=")+2);
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

	/**\param &Coordinates This is a pointer to the 4D vector which has the
	coordinated of all vertices */
	vector<vector<vector<vector<double> > > > Coordinates;
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
	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;
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
	
	// After calling ghostcell function volumes of all the ghost cells will be decided 
	ghostcell(Coordinates,iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,
	CellVolume, i0GhostCellVolume,j0GhostCellVolume, k0GhostCellVolume,
	 iNiGhostCellVolume,jNjGhostCellVolume,kNkGhostCellVolume,
	  Ni, Nj, Nk);

	// After this point the Ni, Nj, Nk has been decided 
	cout << "Ni, Nj, Nk :-> "<< Ni << "  " << Nj << "  " << Nk << endl;

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
	 Ni, Nj, Nk);
		
	ofstream kullu_mass ;
	kullu_mass.open("./Results/outputfiles/Residual.csv");
	// kullu_mass <<  "t(secs)" << "," << "DensityResidual"  << "," <<
	// "xMomentumResidual" << "," << "yMomentumResidual" <<","<< 
	// "zMomentumResidual" << "," << "EnergyResidual" << endl ;

	//Intermidian variables used by the solver
	double iCellInterfaceVolume;
	/**\param iCellInterfaceVolume Average of 
	right and left cell volume in i direction*/
	double jCellInterfaceVolume; 
	/**\param jCellInterfaceVolume Average of 
	right and left cell volume in j direction*/
	double kCellInterfaceVolume; 
	/**\param kCellInterfaceVolume Average of 
	right and left cell volume in k direction*/
	
	double Density ;
	double Pressure ;
	double Velocity ;
	double VelocitySound;

	vector <double> LeftConservedVariables(5) ;
	vector <double> RightConservedVariables(5) ;
	
	/**\param i0GhostConservedVariable Ghost cell Conserved variables array at i = 0*/
	/**\param j0GhostConservedVariable Ghost cell Conserved variables array at i = 0*/
	/**\param k0GhostConservedVariable Ghost cell Conserved variables array at i = 0*/
	/**\param iNiGhostConservedVariable Ghost cell Conserved variables array at i = Ni*/
	/**\param jNjGhostConservedVariable Ghost cell Conserved variables array at j = Nj*/
	/**\param kNkGhostConservedVariable Ghost cell Conserved variables array at k = Nk*/
	Dim4 i0GhostConservedVariable(1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	Dim4 iNiGhostConservedVariable(1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	Dim4 j0GhostConservedVariable(Ni,Dim3(1,Dim2(Nk,Dim1(5))));
	Dim4 jNjGhostConservedVariable(Ni,Dim3(1,Dim2(Nk,Dim1(5))));
	Dim4 k0GhostConservedVariable(Ni,Dim3(Nj,Dim2(1,Dim1(5))));
	Dim4 kNkGhostConservedVariable(Ni,Dim3(Nj,Dim2(1,Dim1(5))));
	
	// Iterations starts here 
	// This file is opened to store the residuals at each time step
	for (int t = 0; t < TotalIteration; ++t)
	{
		if(TimeSteping == "Global")
		{
			deltat = 1000.0;
			// Updating the delta t after every iteration
			for (int i = 0; i < Ni; ++i)
			{
				for (int j = 0; j < Nj; ++j)
				{
					for (int k = 0; k < Nk; ++k)
					{
						/**\bug Local time step needs to be used to reduce the 
						simulation time*/
						Density = ConservedVariables[i][j][k][0];
						Pressure = (SpecificHeatRatio -1)*( ConservedVariables[i][j][k][4] - 0.5*(
						pow(ConservedVariables[i][j][k][1],2)+pow(ConservedVariables[i][j][k][2],2)+
						pow(ConservedVariables[i][j][k][3],2))/Density ) ;  

						Velocity = sqrt(pow(ConservedVariables[i][j][k][1],2)+
							pow(ConservedVariables[i][j][k][2],2)+pow(ConservedVariables[i][j][k][3],2))/Density;
						
						VelocitySound = sqrt(SpecificHeatRatio*Pressure/Density);
						if(deltat > (CFL*delta_s[i][j][k])/(Velocity+VelocitySound))
						{
							deltat = (CFL*delta_s[i][j][k])/(Velocity+VelocitySound);
						}
					}
				}
			}
		}
		
		// cout << "timestep  = " << t << endl ; 
		// Before every time step we need to have proper value in the ghost 
		// cells So, BC takes  care of Inlet, Exit, y-wall and Z-wall BC
		BC(ConservedVariables,
		iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,
		i0GhostConservedVariable,j0GhostConservedVariable,
		k0GhostConservedVariable,iNiGhostConservedVariable,
		jNjGhostConservedVariable,kNkGhostConservedVariable,
		Ni, Nj, Nk);
		

		// i faces calculation
		for (int i = 0; i < Ni+1; ++i)
		{
			for (int j =0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{	
					// i interface volume
					if(i == 0)
					{
						iCellInterfaceVolume = 0.5*(i0GhostCellVolume[0][j][k] + 
							CellVolume[0][j][k]);  
						LeftConservedVariables = i0GhostConservedVariable[0][j][k];
						RightConservedVariables = ConservedVariables[0][j][k];

						// Calculating the flux at the -0.5 interface 
						netfluxAUSM irightface(LeftConservedVariables,
						RightConservedVariables,
						iFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[0][j][k][l] +=(deltat/
								iCellInterfaceVolume)*(irightface.NetFlux[l]);
						}
					}
					#if 1
					else if(i == Ni)
					{
						iCellInterfaceVolume = 0.5*(CellVolume[Ni-1][j][k] +
							iNiGhostCellVolume[0][j][k]);
						LeftConservedVariables = ConservedVariables[Ni-1][j][k];
						RightConservedVariables = iNiGhostConservedVariable[0][j][k];
						
						// Calculating the flux at the Ni-0.5 interface 
						netfluxAUSM irightface(LeftConservedVariables,
						RightConservedVariables,
						iFaceAreaVector[i][j][k]);
						
						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[Ni-1][j][k][l] -=(deltat/
								iCellInterfaceVolume)*(irightface.NetFlux[l]);
						}
						
					}
					else
					{
						iCellInterfaceVolume = 0.5*(CellVolume[i-1][j][k] + 
						CellVolume[i][j][k]);
						LeftConservedVariables = ConservedVariables[i-1][j][k];
						RightConservedVariables = ConservedVariables[i][j][k];

						// Calculating the flux at the i-0.5 interface 
						netfluxAUSM irightface(LeftConservedVariables,
						RightConservedVariables,
						iFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[i-1][j][k][l] -=(deltat/
								iCellInterfaceVolume)*(irightface.NetFlux[l]);
							ConservedVariablesNew[i][j][k][l] +=(deltat/
								iCellInterfaceVolume)*(irightface.NetFlux[l]);
						}

					}
					#endif
				}
			}
		}

		// j faces calculation
		for (int i = 0; i < Ni; ++i)
		{
			for (int j =0; j < Nj+1; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{	
					if(j == 0)
					{
						jCellInterfaceVolume = 0.5*(j0GhostCellVolume[i][0][k] + 
							CellVolume[i][0][k]);
						LeftConservedVariables = j0GhostConservedVariable[i][0][k];
						RightConservedVariables = ConservedVariables[i][0][k];

						// Calculating the flux at the -0.5 interface 
						netfluxAUSM jrightface(LeftConservedVariables,
						RightConservedVariables,
						jFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[i][0][k][l] +=(deltat/
								jCellInterfaceVolume)*(jrightface.NetFlux[l]);
						}
					}
					else if(j == Nj)
					{
						jCellInterfaceVolume = 0.5*(CellVolume[i][Nj][k] +
						jNjGhostCellVolume[i][0][k]);
						LeftConservedVariables = ConservedVariables[i][Nj-1][k];
						RightConservedVariables = jNjGhostConservedVariable[i][0][k];

						// Calculating the flux at the Nj-0.5 interface 
						netfluxAUSM jrightface(LeftConservedVariables,
						RightConservedVariables,
						jFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[i][Nj-1][k][l] -=(deltat/
								jCellInterfaceVolume)*(jrightface.NetFlux[l]);
						}
					}
					else
					{
						jCellInterfaceVolume = 0.5*(CellVolume[i][j-1][k] + 
						CellVolume[i][j][k]);
						LeftConservedVariables = ConservedVariables[i][j-1][k];
						RightConservedVariables = ConservedVariables[i][j][k];

						// Calculating the flux at the j-0.5 interface 
						netfluxAUSM jrightface(LeftConservedVariables,
						RightConservedVariables,
						jFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[i][j-1][k][l] -=(deltat/
								jCellInterfaceVolume)*(jrightface.NetFlux[l]);
							ConservedVariablesNew[i][j][k][l] +=(deltat/
								jCellInterfaceVolume)*(jrightface.NetFlux[l]);
						}
					}

				}
			}
		}

		
		for (int i = 0; i < Ni; ++i)
		{
			for (int j =0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk+1; ++k)
				{
					// k interface volume
					if(k == 0)
					{
						kCellInterfaceVolume = 0.5*(k0GhostCellVolume[i][j][0] + 
							CellVolume[i][j][0]);
						LeftConservedVariables = k0GhostConservedVariable[i][j][0];
						RightConservedVariables = ConservedVariables[i][j][0];

						// Calculating the flux at the -0.5 interface 
						netfluxAUSM krightface(LeftConservedVariables,
						RightConservedVariables,
						kFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[i][j][0][l] +=(deltat/
								kCellInterfaceVolume)*(krightface.NetFlux[l]);
						}
					}
					
					else if(k == Nk)
					{
						kCellInterfaceVolume = 0.5*(CellVolume[i][j][Nk-1]+
							kNkGhostCellVolume[i][j][0]);
						LeftConservedVariables = ConservedVariables[i][j][Nk-1];
						RightConservedVariables = kNkGhostConservedVariable[i][j][0];

						// Calculating the flux at the Nk-0.5 interface 
						netfluxAUSM krightface(LeftConservedVariables,
						RightConservedVariables,
						kFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
						#if 0
							ConservedVariablesNew[i][j][Nk-1][l] -=(deltat/
								kCellInterfaceVolume)*(krightface.NetFlux[l]);
						#endif
						}
					}
					else
					{
						jCellInterfaceVolume = 0.5*(CellVolume[i][j][k-1] + 
						CellVolume[i][j][k]);
						LeftConservedVariables = ConservedVariables[i][j][k-1];
						RightConservedVariables = ConservedVariables[i][j][k];

						// Calculating the flux at the j-0.5 interface 
						netfluxAUSM krightface(LeftConservedVariables,
						RightConservedVariables,
						kFaceAreaVector[i][j][k]);

						for (int l = 0; l < 5; ++l)
						{
							ConservedVariablesNew[i][j][k-1][l] -=(deltat/
								kCellInterfaceVolume)*(krightface.NetFlux[l]);
							ConservedVariablesNew[i][j][k][l] +=(deltat/
								kCellInterfaceVolume)*(krightface.NetFlux[l]);
						}
					}
				}
			}
		}					


		// Residual calculation after each time step and writing the all 
		//residuals into the file
		double DensityResidual = 0.0 ; 
		/**\param DensityResidual Density residual*/
		double xMomentumResidual = 0.0 ; 
		/**\param xMomentumResidual x Momentum residual*/
		double yMomentumResidual = 0.0 ; 
		/**\param yMomentumResidual y Momentum residual*/
		double zMomentumResidual = 0.0 ; 
		/**\param zMomentumResidual z Momentum residual*/
		double EnergyResidual = 0.0 ; 
		/**\param Energy residual */

		int TotalGridPoints = Ni*Nj*Nk ; 
		for (int x = 0; x < Ni; ++x)
		{
			for (int y = 0; y < Nj; ++y)
			{
				DensityResidual   += pow((ConservedVariablesNew[x][y][1][0] -
					ConservedVariables[x][y][1][0]),2);
				xMomentumResidual += pow((ConservedVariablesNew[x][y][1][1] - 
					ConservedVariables[x][y][1][1]),2);     
				yMomentumResidual += pow((ConservedVariablesNew[x][y][1][2] - 
					ConservedVariables[x][y][1][2]),2);     
				zMomentumResidual += pow((ConservedVariablesNew[x][y][1][3] - 
					ConservedVariables[x][y][1][3]),2);     
				EnergyResidual    += pow((ConservedVariablesNew[x][y][1][4] - 
					ConservedVariables[x][y][1][4]),2);    
			}
		}
		if (t%10==0)
		{
			// cout << "TotalGridPoints" << TotalGridPoints << endl ;
			kullu_mass << t << "," << t*deltat << "," << sqrt(DensityResidual/
			(Ni*Nj))  << "," << sqrt(xMomentumResidual/(Ni*Nj))
			<< "," << sqrt(yMomentumResidual/(Ni*Nj)) <<","<< sqrt(
			zMomentumResidual/(Ni*Nj)) << "," << sqrt(EnergyResidual/
			(Ni*Nj)) << endl ;			
			
			/*This is to stop the simulation automatically if nan occurs*/
			if(isnan(sqrt(DensityResidual))==1)
			{
				return 0;
			}
		}

		if (t%10==0)
		{
			cout <<  t << "  --->  " << "  "<< deltat << "  "<<  sqrt(DensityResidual) << endl ;
		}

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

		if (t%10 == 0)
		{
			// storing the all conserved variables in one plane
			ofstream kullu_2D ;
			kullu_2D.open("./Results/outputfiles/2D_parameters_B.csv");
			// kullu_2D << "density" << "," << "density*u" << ","<< "density*v"
			// << "," << "density*w" << "," << "energy"  << endl ;
			for (int i = 0; i < Ni; ++i)
			{
				for (int j = 0; j < Nj; ++j)
				{
					kullu_2D << ConservedVariables[i][j][Nk/2][0] << "," << 
					ConservedVariables[i][j][Nk/2][1] <<","<< 
					ConservedVariables[i][j][Nk/2][2] << "," <<
					ConservedVariables[i][j][Nk/2][3] << "," <<
					ConservedVariables[i][j][Nk/2][4] << endl ;
				}
			}
		}

		#if 0
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
		#endif
	} 
	// time progression ends here 
	time(&EndTime) ;
	double SimulationTotalTime = difftime (EndTime,StartTime);
	cout << "Time taken by the solver in secs = " << SimulationTotalTime<<endl;
	// return 0;
	return 0;
}
