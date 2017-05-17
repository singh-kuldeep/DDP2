#ifndef BC_H
#define BC_H

/*! \file BC.h
    \brief This header file implements all three boundary conditions. 
    - Inlet
    - Exit and 
    - Wall boundary

    \author Kuldeep Singh
    \date 2017 
*/

// #include "iostream"
#include "math.h"
#include <vector>
#include <fstream>
using namespace std;
# define RealGasConstant 287.17 

/** \brief Changes the input vector into the unit normal vector.
*\param areaVector A 3D vector.
*\param vectorMagnitude Magnitude of the 3D vector.
*\return void
*/
void getNormal(vector<double> & UnitNormal, vector<double> areaVector)
{
	double vectorMagnitude = sqrt(areaVector[0]*areaVector[0] + areaVector[1]*
		areaVector[1] + areaVector[2]*areaVector[2]);
	UnitNormal[0] = areaVector[0]/vectorMagnitude;
	UnitNormal[1] = areaVector[1]/vectorMagnitude;
	UnitNormal[2] = areaVector[2]/vectorMagnitude;
}

double getMachfromPressureRatio(double Pressure, double TotalPressure, 
	double SpecificHeatRatio)
{
	// double Mach = 
	double Mach = sqrt((2/(SpecificHeatRatio-1))*
	(pow((Pressure/TotalPressure),(-(SpecificHeatRatio-1)/SpecificHeatRatio))-1));
	// cout << Pressure <<","<< TotalPressure << endl;
	return Mach; 
}
/** \brief This function implements the wall boundary condition
*\param AreaVectors Surface faces area vectors.
*\param LiveCellConservedVariables Conserved variables array for the live cell.
*\param GhostCellConservedVariables Conserved variables array for the ghost cell.
*\return void
*/
void WallBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	vector<double> AreaVectors, double SpecificHeatRatio)
{
	std::vector<double> n(3); // Unit normal vector to the face
	getNormal(n,AreaVectors);
	double densityLive = LiveCellConservedVariables[0]; // density in live cell
	double uLive = LiveCellConservedVariables[1]/densityLive;
	double vLive = LiveCellConservedVariables[2]/densityLive;
	double wLive = LiveCellConservedVariables[3]/densityLive;
	double pressureLive = (SpecificHeatRatio-1)*(LiveCellConservedVariables[4]-
		0.5*densityLive*(uLive*uLive + vLive*vLive + wLive*wLive));

	double densityGhost = densityLive ;
	double ughost = (1-2*n[0]*n[0])*uLive + (-2*n[0]*n[1])*vLive + (-2*n[0]*n[2])*wLive;
	double vghost = (-2*n[1]*n[0])*uLive + (1-2*n[1]*n[1])*vLive + (-2*n[1]*n[2])*wLive;
	double wghost = (-2*n[2]*n[0])*uLive + (-2*n[2]*n[1])*vLive + (1-2*n[2]*n[2])*wLive;
	double PressureGhost = pressureLive;


	GhostCellConservedVariables[0] = densityGhost;
	GhostCellConservedVariables[1] = densityGhost*ughost; 
	GhostCellConservedVariables[2] = densityGhost*vghost; 
	GhostCellConservedVariables[3] = densityGhost*wghost;
	GhostCellConservedVariables[4] = (PressureGhost/(SpecificHeatRatio-1)) + 
	0.5*densityGhost*(ughost*ughost+vghost*vghost+wghost*wghost);  
}

/*Four properties are specified and one is extrapolated, based on analysis of
information propagation along characteristic directions in the calculation 
domain
Here density and velocity are imposed and pressure will be a extrapolated  
*/
void SubSonicInletBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	double InletTotalPressure, double InletTotalTemperature, 
	double SpecificHeatRatio)
{
	double densityLive = LiveCellConservedVariables[0]; // density in live cell
	double uLive = LiveCellConservedVariables[1]/densityLive;
	double vLive = LiveCellConservedVariables[2]/densityLive;
	double wLive = LiveCellConservedVariables[3]/densityLive;
	double pressureLive = (SpecificHeatRatio-1)*(LiveCellConservedVariables[4]-
		0.5*densityLive*(uLive*uLive + vLive*vLive + wLive*wLive));

	double PressureGhost = pressureLive ; // simple extrapolation 
	
	double MachGhost = getMachfromPressureRatio(PressureGhost,
		InletTotalPressure,SpecificHeatRatio);

	double TemperatureGhost = InletTotalTemperature*
	pow((1+((SpecificHeatRatio-1)*pow(MachGhost,2)/2)),-1.0); 
	
	double Density = PressureGhost/(RealGasConstant*TemperatureGhost);

	double SoundSpeed = sqrt(SpecificHeatRatio*RealGasConstant*TemperatureGhost);
	double XVelocity = MachGhost*SoundSpeed;
	double YVelocity = 0.0;
	double ZVelocity = 0.0;

	GhostCellConservedVariables[0] = Density;
	GhostCellConservedVariables[1] = Density*XVelocity;
	GhostCellConservedVariables[2] = Density*YVelocity;
	GhostCellConservedVariables[3] = Density*ZVelocity;
	GhostCellConservedVariables[4] = (PressureGhost/(SpecificHeatRatio-1)) + 
	0.5*Density*(XVelocity*XVelocity+YVelocity*YVelocity+ZVelocity*ZVelocity);
}

// Only one quantity is imposed 
void SubSonicExitBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	double ExitPressure, double SpecificHeatRatio)
{
	double densityLive = LiveCellConservedVariables[0]; // density in live cell
	double uLive = LiveCellConservedVariables[1]/densityLive;
	double vLive = LiveCellConservedVariables[2]/densityLive;
	double wLive = LiveCellConservedVariables[3]/densityLive;
	
	// Four parameters are extrapolated form lie cell
	double densityGhost = densityLive ;
	double ughost = uLive;
	double vghost = vLive;
	double wghost = wLive;

	// Exit Pressure is imposed 
	double PressureGhost = ExitPressure; 

	GhostCellConservedVariables[0] = densityGhost;
	GhostCellConservedVariables[1] = densityGhost*ughost; 
	GhostCellConservedVariables[2] = densityGhost*vghost; 
	GhostCellConservedVariables[3] = densityGhost*wghost;
	GhostCellConservedVariables[4] = (PressureGhost/(SpecificHeatRatio-1)) + 
	0.5*densityGhost*(ughost*ughost+vghost*vghost+wghost*wghost);	
}

// all the parameters are extrapolated 
void SuperSonicExitBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables)
{
	for (int l = 0; l < 5; ++l)
	 {
	 	GhostCellConservedVariables[l] = LiveCellConservedVariables[l];
	 } 
}

void SuperSonicInletBC(vector<double> & GhostCellConservedVariables,
	double InletTotalPressure, double InletTotalTemperature, double InletMach,
	double SpecificHeatRatio)
{
	double MachGhost = InletMach; 

	double PressureGhost = InletTotalPressure*
	pow((1+((SpecificHeatRatio-1)*pow(MachGhost,2)/2)),(-SpecificHeatRatio/(SpecificHeatRatio-1)));  
	

	double TemperatureGhost = InletTotalTemperature*
	pow((1+((SpecificHeatRatio-1)*pow(MachGhost,2)/2)),-1.0); 
	
	double Density = PressureGhost/(RealGasConstant*TemperatureGhost);

	double SoundSpeed = sqrt(SpecificHeatRatio*RealGasConstant*TemperatureGhost);

	double XVelocity = MachGhost*SoundSpeed;
	double YVelocity = 0.0;
	double ZVelocity = 0.0;

	GhostCellConservedVariables[0] = Density;
	GhostCellConservedVariables[1] = Density*XVelocity;
	GhostCellConservedVariables[2] = Density*YVelocity;
	GhostCellConservedVariables[3] = Density*ZVelocity;
	GhostCellConservedVariables[4] = (PressureGhost/(SpecificHeatRatio-1)) + 
	0.5*Density*(XVelocity*XVelocity+YVelocity*YVelocity+ZVelocity*ZVelocity);
}

/** \brief Function BC() implements the boundary condition. 
* Here two ghost cell are used to implement the boundary condition. In simple 
words this function calculates the conserved variables for all ghost cells.
*For inlet it uses the stagnation parameters, for exit it simply uses the live 
cell parameters and copies them into the ghost cells, and for wall boundary
* it uses the fact that flow should be parallel to the wall. 
*\param [in] ConservedVariables This is the pointer to the 4D vector where all 
the conserved variables of previous time step are stored.
*\param [in] &iFaceAreaVector This is a pointer to the 4D vector which has the 
area vector of all faces which are in "i" direction.  
*\param [in] &jFaceAreaVector This is a pointer to the 4D vector which has the 
area vector of all faces which are in "j" direction.  
*\param [in] &kFaceAreaVector This is a pointer to the 4D vector which has the 
area vector of all faces which are in "k" direction.  
*\param [in] Ni Number of cells in in "i" direction.  
*\param [in] Nj Number of cells in in "j" direction.  
*\param [in] Nk Number of cells in in "k" direction.  
*\return void
*/

void BC(
	vector<vector<vector<vector<double> > > > ConservedVariables,
	vector<vector<vector<vector<double> > > > iFaceAreaVector,
	vector<vector<vector<vector<double> > > > jFaceAreaVector,
	vector<vector<vector<vector<double> > > > kFaceAreaVector,
	
	vector<vector<vector<vector<double> > > > & i0GhostConservedVariable,
	vector<vector<vector<vector<double> > > > & j0GhostConservedVariable,
	vector<vector<vector<vector<double> > > > & k0GhostConservedVariable,
	vector<vector<vector<vector<double> > > > & iNiGhostConservedVariable,
	vector<vector<vector<vector<double> > > > & jNjGhostConservedVariable,
	vector<vector<vector<vector<double> > > > & kNkGhostConservedVariable,
	int Ni, int Nj, int Nk, double SpecificHeatRatio)
{	
	// if premetive variables are given at the inlet
	double InletDensity ; 
	double InletXVelocity ; 
	double InletYVelocity ; 
	double InletZVelocity ; 
	double InletStaticPressure ;

	double ExitStaticPressure ;

	// if total quantities are given at the inlet 
	double InletTotalTemperature;	
	double InletTotalPressure;

	double InletMach; //only in case of supersonic inlet

	string BoundaryConditionati0 ;
	string BoundaryConditionatj0 ;
	string BoundaryConditionatk0 ;
	string BoundaryConditionatiNi ;
	string BoundaryConditionatjNj ;
	string BoundaryConditionatkNk ;

	ifstream infile("inputfile");
	string aline;

	// reading the input file 
	while(!infile.eof())// file ended
	{
		getline(infile,aline); // reading line form file

		if (aline.find( "//" )!=0 && aline.empty()==false) 
		{
			if(aline.find("InletDensity")!=string::npos)
			{
				InletDensity = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletXVelocity")!=string::npos)
			{
				InletXVelocity = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletYVelocity")!=string::npos)
			{
				InletYVelocity = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletZVelocity")!=string::npos)
			{
				InletZVelocity = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletStaticPressure")!=string::npos)
			{
				InletStaticPressure = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("ExitStaticPressure")!=string::npos)
			{
				ExitStaticPressure = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletTotalTemperature")!=string::npos)
			{
				InletTotalTemperature = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletTotalPressure")!=string::npos)
			{
				InletTotalPressure = stod (aline.substr(aline.find("=")+1));
			}
			else if(aline.find("InletMach")!=string::npos)
			{
				InletMach = stod (aline.substr(aline.find("=")+1));
			}
			else if (aline.find("BoundaryConditionati0")!=string::npos)
			{
				BoundaryConditionati0 = aline.substr(aline.find("=")+2); 
			}
			else if (aline.find("BoundaryConditionatj0")!=string::npos)
			{
				BoundaryConditionatj0 = aline.substr(aline.find("=")+2); 
			}
			else if (aline.find("BoundaryConditionatk0")!=string::npos)
			{
				BoundaryConditionatk0 = aline.substr(aline.find("=")+2); 
			}
			else if (aline.find("BoundaryConditionatiNi")!=string::npos)
			{
				BoundaryConditionatiNi = aline.substr(aline.find("=")+2); 
			}
			else if (aline.find("BoundaryConditionatjNj")!=string::npos)
			{
				BoundaryConditionatjNj = aline.substr(aline.find("=")+2); 
			}
			else if (aline.find("BoundaryConditionatkNk")!=string::npos)
			{
				BoundaryConditionatkNk = aline.substr(aline.find("=")+2); 
			}				
		}
	}

	/* Implementing the boundary condition based on the options given 
	in input file */
	for (int j = 0; j < Nj; ++j)
	{
		for (int k = 0; k < Nk; ++k)
		{
			/*i0*/
			if(BoundaryConditionati0 == "Wall")
			{
				WallBC(i0GhostConservedVariable[0][j][k],
				ConservedVariables[0][j][k],iFaceAreaVector[0][j][k],
				SpecificHeatRatio);
			}
			else if(BoundaryConditionati0 == "SubSonicInlet")
			{
				SubSonicInletBC(i0GhostConservedVariable[0][j][k],
				ConservedVariables[0][j][k],
				InletTotalPressure, InletTotalTemperature,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionati0 == "SuperSonicInlet")
			{
				SuperSonicInletBC(i0GhostConservedVariable[0][j][k],
				InletTotalPressure, InletTotalTemperature, InletMach,
				SpecificHeatRatio);	
			}
			else if(BoundaryConditionati0 == "SubSonicExit")
			{
				SubSonicExitBC(i0GhostConservedVariable[0][j][k],
				ConservedVariables[0][j][k], ExitStaticPressure,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionati0 == "SuperSonicExit")
			{
				SuperSonicExitBC(i0GhostConservedVariable[0][j][k],
				ConservedVariables[0][j][k]);
			}

			/*iNi*/
			if(BoundaryConditionatiNi == "Wall")
			{
				WallBC(iNiGhostConservedVariable[0][j][k],
				ConservedVariables[Ni-1][j][k],iFaceAreaVector[Ni][j][k],
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatiNi == "SubSonicInlet")
			{
				SubSonicInletBC(iNiGhostConservedVariable[0][j][k],
				ConservedVariables[Ni-1][j][k],
				InletTotalPressure, InletTotalTemperature,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatiNi == "SuperSonicInlet")
			{
				SuperSonicInletBC(iNiGhostConservedVariable[0][j][k],
				InletTotalPressure, InletTotalTemperature,InletMach,
				SpecificHeatRatio);	
			}
			else if(BoundaryConditionatiNi == "SubSonicExit")
			{
				SubSonicExitBC(iNiGhostConservedVariable[0][j][k],
				ConservedVariables[Ni-1][j][k],
				ExitStaticPressure,SpecificHeatRatio);
			}
			else if(BoundaryConditionatiNi == "SuperSonicExit")
			{
				SuperSonicExitBC(iNiGhostConservedVariable[0][j][k],
				ConservedVariables[Ni-1][j][k]);
			}
		}
	}

	for (int i = 0; i < Ni; ++i)
	{
		for (int k = 0; k < Nk; ++k)
		{
			/*j0*/
			if(BoundaryConditionatj0 == "Wall")
			{
				WallBC(j0GhostConservedVariable[i][0][k],
				ConservedVariables[i][0][k],jFaceAreaVector[i][0][k],
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatj0 == "SubSonicInlet")
			{
				SubSonicInletBC(j0GhostConservedVariable[i][0][k],
				ConservedVariables[i][0][k],
				InletTotalPressure, InletTotalTemperature,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatj0 == "SuperSonicInlet")
			{
				SuperSonicInletBC(j0GhostConservedVariable[i][0][k],
				InletTotalPressure, InletTotalTemperature,InletMach,		
				SpecificHeatRatio);	
			}
			else if(BoundaryConditionatj0 == "SubSonicExit")
			{
				SubSonicExitBC(j0GhostConservedVariable[i][0][k],
				ConservedVariables[i][0][k],
				ExitStaticPressure,SpecificHeatRatio);
			}
			else if(BoundaryConditionatj0 == "SuperSonicExit")
			{
				SuperSonicExitBC(j0GhostConservedVariable[i][0][k],
				ConservedVariables[i][0][k]);
			}

			/*jNj*/
			if(BoundaryConditionatjNj == "Wall")
			{
				WallBC(jNjGhostConservedVariable[i][0][k],
				ConservedVariables[i][Nj-1][k],jFaceAreaVector[i][Nj][k],
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatjNj == "SubSonicInlet")
			{
				SubSonicInletBC(jNjGhostConservedVariable[i][0][k],
				ConservedVariables[i][Nj-1][k],
				InletTotalPressure, InletTotalTemperature,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatjNj == "SuperSonicInlet")
			{
				SuperSonicInletBC(jNjGhostConservedVariable[i][0][k],
				InletTotalPressure,InletTotalTemperature,InletMach,
				SpecificHeatRatio);	
			}
			else if(BoundaryConditionatjNj == "SubSonicExit")
			{
				SubSonicExitBC(jNjGhostConservedVariable[i][0][k],
				ConservedVariables[i][Nj-1][k],
				ExitStaticPressure,SpecificHeatRatio);
			}
			else if(BoundaryConditionatjNj == "SuperSonicExit")
			{
				SuperSonicExitBC(jNjGhostConservedVariable[i][0][k],
				ConservedVariables[i][Nj-1][k]);
			}
		}
	}

	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			/*k0*/
			if(BoundaryConditionatk0 == "Wall")
			{
				WallBC(k0GhostConservedVariable[i][j][0],
				ConservedVariables[i][j][0],kFaceAreaVector[i][j][0],
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatk0 == "SubSonicInlet")
			{
				SubSonicInletBC(k0GhostConservedVariable[i][j][0],
				ConservedVariables[i][j][0],
				InletTotalPressure, InletTotalTemperature,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatk0 == "SuperSonicInlet")
			{
				SuperSonicInletBC(k0GhostConservedVariable[i][j][0],
				InletTotalPressure, InletTotalTemperature,InletMach,
				SpecificHeatRatio);	
			}
			else if(BoundaryConditionatk0 == "SubSonicExit")
			{
				SubSonicExitBC(k0GhostConservedVariable[i][j][0],
				ConservedVariables[i][j][0],
				ExitStaticPressure,SpecificHeatRatio);
			}
			else if(BoundaryConditionatk0 == "SuperSonicExit")
			{
				SuperSonicExitBC(k0GhostConservedVariable[i][j][0],
				ConservedVariables[i][j][0]);
			}

			/*kNk*/
			if(BoundaryConditionatkNk == "Wall")
			{
				WallBC(kNkGhostConservedVariable[i][j][0],
				ConservedVariables[i][j][Nk-1],kFaceAreaVector[i][j][Nk],
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatkNk == "SubSonicInlet")
			{
				SubSonicInletBC(kNkGhostConservedVariable[i][j][0],
				ConservedVariables[i][j][Nk-1],
				InletTotalPressure, InletTotalTemperature,
				SpecificHeatRatio);
			}
			else if(BoundaryConditionatkNk == "SuperSonicInlet")
			{
				SuperSonicInletBC(kNkGhostConservedVariable[i][j][0],
				InletTotalPressure, InletTotalTemperature,InletMach,
				SpecificHeatRatio);	
			}
			else if(BoundaryConditionatkNk == "SubSonicExit")
			{
				SubSonicExitBC(kNkGhostConservedVariable[i][j][0],
				ConservedVariables[i][j][Nk-1],
				ExitStaticPressure,SpecificHeatRatio);
			}
			else if(BoundaryConditionatkNk == "SuperSonicExit")
			{
				SuperSonicExitBC(kNkGhostConservedVariable[i][j][0],
				ConservedVariables[i][j][Nk-1]);
			}
		}
	}

}
#endif // BC.h ends here
