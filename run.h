/*! \file run.h
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
#include "grid_nozzle.h"
#include "BC.h"
// #include "dt.h"
// #include "grid_straight_duct.h"
// #include "grid_bump.h"
// #include "grid_diverging_duct.h"

using namespace std ;


/** @brief This function implements the boundary condition, iFaceAreaVector is not required Because currently the flow in x direction and 2D flow */
void BC(vector<vector<vector<vector<double> > > > & ConservedVariables, vector<vector<vector<vector<double> > > > &
 jFaceAreaVector, vector<vector<vector<vector<double> > > > & kFaceAreaVector, int Ni, int Nj, int Nk) ;

/** @brief This function generates the area vector and cell volumes inside the domain whole domain*/
void grid(vector<vector<vector<vector<double> > > > & iFaceAreaVector , 
	vector<vector<vector<vector<double> > > > & jFaceAreaVector,
	vector<vector<vector<vector<double> > > > & kFaceAreaVector, 
	vector<vector<vector<double> > >& CellVolume,
	vector<vector<vector<double> > >& delta_s,
	int & Ni, int & Nj, int & Nk);

/*! \brief This function runs the solver.
    \warning Currently not using this, because  grid() is not calculating ds value properly. 
    So recheck this function as well after fixing the grid() function.
	\return double
*/
void run ()
{
	time_t StartTime; /**\param StartTime Simulation starting time*/
	time_t EndTime ; /**\param EndTime Simulation ending time*/
	time(&StartTime); // noting the starting time

	double DeltaT = 0.00000015; /**\param DeltaT Time step*/
	double TIME = 10000000*DeltaT;
	int IterationValues = 100000 /**\param IterationValues Total iterations = floor(TIME/DeltaT)*/ ;

	int Ni;/**\param Ni Number of cells in in "i" direction.*/
	int Nj;/**\param Nj Number of cells in in "j" direction.*/
	int Nk;/**\param Nk Number of cells in in "k" direction.*/
	
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
	int Nk = 1+4; // Because this is 2D-simulation so no need to take large number of grids in z direction 
	#endif  

	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> matrix4D;

	/**\param &iFaceAreaVector This is a pointer to the 4D vector which has the area vector of all faces which are in "i" direction.*/
	std::vector<std::vector<std::vector<std::vector<double> > > > iFaceAreaVector ;
	/**\param &jFaceAreaVector This is a pointer to the 4D vector which has the area vector of all faces which are in "j" direction.*/
	std::vector<std::vector<std::vector<std::vector<double> > > > jFaceAreaVector ;
	/**\param &kFaceAreaVector This is a pointer to the 4D vector which has the area vector of all faces which are in "k" direction.*/
	std::vector<std::vector<std::vector<std::vector<double> > > > kFaceAreaVector ;
	/**\param CellVolumeIn Input pointer to cell volumes*/
	std::vector<std::vector<std::vector<double> > > CellVolume ;
	/**\param delta_s Minimum distance*/
	std::vector<std::vector<std::vector<double> > > delta_s ;
	
	grid(iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,CellVolume,delta_s,Ni,Nj,Nk);

	cout << "Ni, Nj, Nk :-> "<< Ni << "  " << Nj << "  " << Nk << endl;

	// this store previous values of variables 
	matrix4D ConservedVariables(Ni,Dim3(Nj,Dim2(Nk,Dim1(5)))); 
	/**\param ConservedVariables This is the pointer to the 4D vector where all the 
	conserved variables ([Density , x-momentum, y-momentum, z-momentum, Energy]) of previous time step are stored.*/

	/**\param ConservedVariablesNew This is the pointer to the 4D vector where all the 
	conserved variables ([Density , x-momentum, y-momentum, z-momentum, Energy]) of current/new time step are stored.*/
	matrix4D ConservedVariablesNew(Ni,Dim3(Nj,Dim2(Nk,Dim1(5)))); 
	

	// Initial conditions(these are just random values ) to start 
	for (int i =0; i < Ni; ++i)
	{
		for (int j =0; j < Nj; ++j)
		{
			for (int k =0; k < Nk; ++k)
			{
				ConservedVariables[i][j][k][0] = 1.16;
				ConservedVariables[i][j][k][1] = 0 ;
				ConservedVariables[i][j][k][2] = 0 ;
				ConservedVariables[i][j][k][3] = 0 ;
				ConservedVariables[i][j][k][4] = 271767;

				ConservedVariablesNew[i][j][k][0] = 1.16;
				ConservedVariablesNew[i][j][k][1] = 0 ;
				ConservedVariablesNew[i][j][k][2] = 0 ;
				ConservedVariablesNew[i][j][k][3] = 0 ;
				ConservedVariablesNew[i][j][k][4] = 271767;
			}
		}
	}

			
	/**@bug Every time simulation starts from first iteration. So, to save the simulation it is good to start from the 
	last solution as the initial condition*/
	#if 0
	cout << "Do you wants to start the simulation, where you left?(enter y)" << endl 
	<< "Other wise press any key"<< endl;
	char oldstart;
	cin>> oldstart;
	if (oldstart=='y')
	{
		ifstream nozzleData ("restart.csv");
		string aline;

		// Initial condition are taken from the previously solved part
		for (int i = 0; i < Ni; ++i)
		{
			for (int j = 0; j < Nj; ++j)
			{
				for (int k = 0; k < Nk; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						getline(nozzleData,aline);
						ConservedVariables[i][j][k][l] = atof(aline.c_str());
						ConservedVariablesNew[i][j][k][l] = ConservedVariables[i][j][k][l];
					}
				}
			}
		}	
	}
	#endif

	ofstream kullu_mass ;
	kullu_mass.open("Residual_Nozzle.csv");
	// kullu_mass <<  "t(secs)" << "," << "DensityResidual"  << "," << "xMomentumResidual" << "," <<
	// 	 "yMomentumResidual" <<","<< "zMomentumResidual" << "," << "EnergyResidual" << endl ;

	double iCellInterfaceVolume; /**\param iCellInterfaceVolume Average of right and left cell volume in i direction*/
	double jCellInterfaceVolume; /**\param jCellInterfaceVolume Average of right and left cell volume in j direction*/
	double kCellInterfaceVolume; /**\param kCellInterfaceVolume Average of right and left cell volume in k direction*/
	
	// Iterations starts here 
	// This file is opened to store the residuals at each time step
	for (int t = 0; t < IterationValues; ++t)
	{
		// cout << "timestep  = " << t << endl ; 
		
		// Before every time step we need to have proper value in the ghost cells 
		// So Boundary condition has been applied by calling the function BC()
		// BC takes  care of Inlet, Exit, y-wall and Z-wall boundary condition
		BC(ConservedVariables,jFaceAreaVector,kFaceAreaVector,Ni,Nj,Nk); 
		
		// next time step ConservedVariabless calculation
		for (int i = 1; i < Ni-2; ++i)
		// Or (int i = 2; i < (Ni+1)-2; ++i) Total Ni+1 interface, 2 used in BC implementation
		{
			for (int j = 1; j < Nj-2; ++j)
			{
				for (int k = 1; k < Nk-2; ++k)
				{	
					/**\bug Local time step needs to be used to reduce the simulation time*/

					// DeltaT = dt(i+1,j+1,2,delta_s,ConservedVariables);
					// cout << "dt   " << DeltaT << endl;

					iCellInterfaceVolume = 0.5*(CellVolume[i][j][k] + CellVolume[i+1][j][k]);
					jCellInterfaceVolume = 0.5*(CellVolume[i][j][k] + CellVolume[i][j+1][k]);
					kCellInterfaceVolume = 0.5*( CellVolume[i][j][k] + CellVolume[i][j][k+1]);

					// net flux using the class netfluxinterface
					netfluxinterface xrightface(ConservedVariables[i-1][j][k],ConservedVariables[i][j][k],ConservedVariables[i+1][j][k],
							ConservedVariables[i+2][j][k] ,iFaceAreaVector[i-1][j][k] ,iFaceAreaVector[i][j][k] ,iFaceAreaVector[i+1][j][k],
						CellVolume[i-1][j][k], CellVolume[i][j][k], CellVolume[i+1][j][k], CellVolume[i+2][j][k] ,DeltaT) ;

					netfluxinterface yrightface(ConservedVariables[i][j-1][k],ConservedVariables[i][j][k],ConservedVariables[i][j+1][k],
							ConservedVariables[i][j+2][k] ,jFaceAreaVector[i][j-1][k] ,jFaceAreaVector[i][j][k] ,jFaceAreaVector[i][j+1][k],
						CellVolume[i][j-1][k], CellVolume[i][j][k], CellVolume[i][j+1][k], CellVolume[i][j+2][k] ,DeltaT) ;

					netfluxinterface zrightface(ConservedVariables[i][j][k-1],ConservedVariables[i][j][k],ConservedVariables[i][j][k+1],
							ConservedVariables[i][j][k+2] ,kFaceAreaVector[i][j][k-1] ,kFaceAreaVector[i][j][k] ,kFaceAreaVector[i][j][k+1],
						CellVolume[i][j][k-1], CellVolume[i][j][k], CellVolume[i][j][k+1], CellVolume[i][j][k+2] ,DeltaT) ;

					// updating the ConservedVariablesNew using flux at the right interfaces
					for (int l = 0; l < 5; ++l)
					{
						ConservedVariablesNew[i][j][k][l] -=(DeltaT/iCellInterfaceVolume)*(xrightface.NetFlux[l]);
						ConservedVariablesNew[i+1][j][k][l] +=(DeltaT/iCellInterfaceVolume)*(xrightface.NetFlux[l]);

						ConservedVariablesNew[i][j][k][l] -=(DeltaT/jCellInterfaceVolume)*(yrightface.NetFlux[l]);
						ConservedVariablesNew[i][j+1][k][l] +=(DeltaT/jCellInterfaceVolume)*(yrightface.NetFlux[l]);

						ConservedVariablesNew[i][j][k][l] -=(DeltaT/kCellInterfaceVolume)*(zrightface.NetFlux[l]);
						ConservedVariablesNew[i][j][k+1][l] +=(DeltaT/kCellInterfaceVolume)*(zrightface.NetFlux[l]);
					}
				}
			}
			
		}

			// Residual calculation after each time step and writing the all residuals into the file
			double DensityResidual = 0.0 ; /**\param DensityResidual Density residual */
			double xMomentumResidual = 0.0 ; /**\param xMomentumResidual x Momentum residual*/
			double yMomentumResidual = 0.0 ; /**\param yMomentumResidual y Momentum residual*/
			double zMomentumResidual = 0.0 ; /**\param zMomentumResidual z Momentum residual*/
			double EnergyResidual = 0.0 ; /**\param Energy residual */

			int TotalGridPoints = 0 ; 
			for (int x = 2; x < Ni-2; ++x)
				{
					for (int y = 2; y < Nj-2; ++y)
					{
						TotalGridPoints   += 1 ; 
						DensityResidual    += pow((ConservedVariablesNew[x][y][2][0] - ConservedVariables[x][y][2][0]),2);
						xMomentumResidual += pow((ConservedVariablesNew[x][y][2][1] - ConservedVariables[x][y][2][1]),2);     
						yMomentumResidual += pow((ConservedVariablesNew[x][y][2][2] - ConservedVariables[x][y][2][2]),2);     
						zMomentumResidual += pow((ConservedVariablesNew[x][y][2][3] - ConservedVariables[x][y][2][3]),2);     
						EnergyResidual     += pow((ConservedVariablesNew[x][y][2][4] - ConservedVariables[x][y][2][4]),2);     
					}
				}
		if (t%10 == 0)
		{
			// cout << "TotalGridPoints" << TotalGridPoints << endl ;
			kullu_mass << t << "," << t*DeltaT << "," << sqrt(DensityResidual/((Ni-4)*(Nj-4)))  << "," << sqrt(xMomentumResidual/((Ni-4)*(Nj-4))) << "," <<
			 sqrt(yMomentumResidual/((Ni-4)*(Nj-4))) <<","<< sqrt(zMomentumResidual/((Ni-4)*(Nj-4))) << "," << sqrt(EnergyResidual/((Ni-4)*(Nj-4))) << endl ;			
		}

		if (t%10==0)
		{
			cout <<  t << "  --->  " << "  "<<  sqrt(DensityResidual) << endl ;
		}

		// before going to the new time step updateing the old conserved variables by new ones.
		for (int i = 2; i < Ni-2; ++i)
		{
			for (int j = 2; j < Nj-2; ++j)
			{
				for (int k = 2; k < Nk-2; ++k)
				{
					for (int l = 0; l < 5; ++l)
					{
						ConservedVariables[i][j][k][l] = ConservedVariablesNew[i][j][k][l] ;
					}				
				}
			}
		}

		if (t%100 == 0)
		{
			// storing the all conserved variables in one plane
			ofstream kullu_2D ;
			kullu_2D.open("2D_parameters_B.csv");
			// kullu_2D << "density" << "," << "density*u" << ","<< "density*v" << "," << "density*w" << "," << "energy"  << endl ;
			for (int i = 2; i < Ni-2; ++i)
			{
				for (int j = 2; j < Nj-2; ++j)
				{
					kullu_2D << ConservedVariables[i][j][Nk/2][0] << "," << ConservedVariables[i][j][Nk/2][1] <<","<< 
					ConservedVariables[i][j][Nk/2][2] << "," << ConservedVariables[i][j][Nk/2][3] << "," << ConservedVariables[i][j][Nk/2][4] << endl ;
				}
			}
		}

		// Restart file is being written to restart the simulation where it was left earlier
		if (t%500 == 0)
		{
			// storing the conserved variables throughout the domain
			ofstream kullu_restart ;
			kullu_restart.open("restart.csv");
			// kullu_restart << "density" << "," << "density*u" << ","<< "density*v" << "," << "density*w" << "," << "energy"  << endl ;
			for (int i = 0; i < Ni; ++i)
			{
				for (int j  = 0; j < Nj; ++j)
				{
					for (int k = 0; k < Nk; ++k)
					{
						kullu_restart << ConservedVariables[i][j][k][0] << endl << ConservedVariables[i][j][k][1] <<endl<< 
						ConservedVariables[i][j][k][2] << endl << ConservedVariables[i][j][k][3] << endl << ConservedVariables[i][j][k][4] << endl ;
					}
				}
			}
		}
	} 
	// time progression ends here 
	time(&EndTime) ;
	double SimulationTotalTime = difftime (EndTime,StartTime);
	cout << "Time taken by the solver in secs = " << SimulationTotalTime << endl ;
	// return 0;
} /**\param */
 
