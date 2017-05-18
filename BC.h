/*! \file BC.h
    \brief Implements the all three boundary conditions. 
    - Inlet
    - Exit and 
    - Wall boundary

    \date 18-May-2017 
*/
#ifndef BC_H
#define BC_H

#include "iostream"
#include "math.h"
#include <vector>
#include <fstream>
#include <cstdlib>

using namespace std;
# define RealGasConstant 287.17 

/*! \fn void getNormal(vector<double> & UnitNormal, vector<double> areaVector)
	\brief Give the unit area vector, using the area vector of the face   
	\param [in] &UnintNormal Pinter of the normal vector
	\param [in] areaVector Area vector of the face  
	\brief Changes the input vector into the unit normal vector.
*/
void getNormal(vector<double> & UnitNormal, vector<double> areaVector)
{
	double vectorMagnitude = sqrt(areaVector[0]*areaVector[0] + areaVector[1]*
		areaVector[1] + areaVector[2]*areaVector[2]);
	UnitNormal[0] = areaVector[0]/vectorMagnitude;
	UnitNormal[1] = areaVector[1]/vectorMagnitude;
	UnitNormal[2] = areaVector[2]/vectorMagnitude;
}

/*! \fn double getMachfromPressureRatio(double Pressure, double TotalPressure, 
	double SpecificHeatRatio)	
	\brief Calculates the Mach number using the total and static pressure
	\param [in] Pressure Static pressure
	\param [in] TotalPressure Total pressure  
	\brief Changes the input vector into the unit normal vector.
	\return Mach number
*/
double getMachfromPressureRatio(double Pressure, double TotalPressure, 
	double SpecificHeatRatio)
{
	double Mach = sqrt((2/(SpecificHeatRatio-1))*
	(pow((Pressure/TotalPressure),(-(SpecificHeatRatio-1)/
		SpecificHeatRatio))-1));
	// cout << Pressure <<","<< TotalPressure << endl;
	return Mach; 
}

/*! \fn  void WallBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	vector<double> AreaVectors, double SpecificHeatRatio)
	\brief This function implements the wall boundary condition
	\param [in] AreaVectors Face area vectors.
	\param [in] LiveCellConservedVariables Conserved variables array for the 
	live cell.
	\param [in,out] GhostCellConservedVariables Conserved variables array for the 
	ghost cell.
*/
void WallBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	vector<double> AreaVectors, double SpecificHeatRatio)
{
	std::vector<double> n(3); // Unit normal vector to the face
	getNormal(n,AreaVectors);
	double densityLive = LiveCellConservedVariables[0];
	double uLive = LiveCellConservedVariables[1]/densityLive;
	double vLive = LiveCellConservedVariables[2]/densityLive;
	double wLive = LiveCellConservedVariables[3]/densityLive;
	double pressureLive = (SpecificHeatRatio-1)*(LiveCellConservedVariables[4]-
		0.5*densityLive*(uLive*uLive + vLive*vLive + wLive*wLive));

	double densityGhost = densityLive ;
	double ughost = (1-2*n[0]*n[0])*uLive + (-2*n[0]*n[1])*vLive + 
	(-2*n[0]*n[2])*wLive;
	double vghost = (-2*n[1]*n[0])*uLive + (1-2*n[1]*n[1])*vLive +
	(-2*n[1]*n[2])*wLive;
	double wghost = (-2*n[2]*n[0])*uLive + (-2*n[2]*n[1])*vLive +
	(1-2*n[2]*n[2])*wLive;
	double PressureGhost = pressureLive;


	GhostCellConservedVariables[0] = densityGhost;
	GhostCellConservedVariables[1] = densityGhost*ughost; 
	GhostCellConservedVariables[2] = densityGhost*vghost; 
	GhostCellConservedVariables[3] = densityGhost*wghost;
	GhostCellConservedVariables[4] = (PressureGhost/(SpecificHeatRatio-1)) + 
	0.5*densityGhost*(ughost*ughost+vghost*vghost+wghost*wghost);  
}

/** \fn void SubSonicInletBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	double InletTotalPressure, double InletTotalTemperature, 
	double SpecificHeatRatio)

	\brief Implements the subsonic inlet boundary condition
	\param [in] LiveCellConservedVariables Conserved variables array for the 
	live cell.
	\param [in,out] GhostCellConservedVariables Conserved variables array for the 
	ghost cell.
	\param [in] InletTotalPressure Total pressure at inlet.
	\param [in] InletTotalTemperature Total temperature at inlet.
	\param [in] SpecificHeatRatio Specific heat ratio 
*/
void SubSonicInletBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	double InletTotalPressure, double InletTotalTemperature, 
	double SpecificHeatRatio)
{
	double densityLive = LiveCellConservedVariables[0];
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

	double SoundSpeed = sqrt(SpecificHeatRatio*RealGasConstant*
		TemperatureGhost);
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

/** \fn void SubSonicExitBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	double ExitPressure, double SpecificHeatRatio)

	\brief Implements the subsonic exit boundary condition
	\param [in] LiveCellConservedVariables Conserved variables array for the 
	live cell.
	\param [in,out] GhostCellConservedVariables Conserved variables array for the 
	ghost cell.
	\param [in] ExitPressure Pressure at exit.
	\param [in] SpecificHeatRatio Specific heat ratio 
*/
void SubSonicExitBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables,
	double ExitPressure, double SpecificHeatRatio)
{
	double densityLive = LiveCellConservedVariables[0];
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

/** \fn void SuperSonicExitBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables)

	\brief Implements the supersonic exit boundary condition, by simple 
	extrapolation
	\param [in] LiveCellConservedVariables Conserved variables array for the 
	live cell.
	\param [in,out] GhostCellConservedVariables Conserved variables array for the 
	ghost cell.
*/
void SuperSonicExitBC(vector<double> & GhostCellConservedVariables,
	vector<double> LiveCellConservedVariables)
{
	for (int l = 0; l < 5; ++l)
	 {
	 	GhostCellConservedVariables[l] = LiveCellConservedVariables[l];
	 } 
}

/** \fn void SuperSonicInletBC(vector<double> & GhostCellConservedVariables,
	double InletTotalPressure, double InletTotalTemperature, double InletMach,
	double SpecificHeatRatio)

	\brief Implements the subsonic exit boundary condition
	\param [in] LiveCellConservedVariables Conserved variables array for the 
	live cell.
	\param [in,out] GhostCellConservedVariables Conserved variables array for the 
	ghost cell.
	\param [in] InletMach Mach number at the inlet.
	\param [in] SpecificHeatRatio Specific heat ratio 
*/
void SuperSonicInletBC(vector<double> & GhostCellConservedVariables,
	double InletTotalPressure, double InletTotalTemperature, double InletMach,
	double SpecificHeatRatio)
{
	double MachGhost = InletMach; 

	double PressureGhost = InletTotalPressure*
	pow((1+((SpecificHeatRatio-1)*pow(MachGhost,2)/2)),(-SpecificHeatRatio/
		(SpecificHeatRatio-1)));  
	

	double TemperatureGhost = InletTotalTemperature*
	pow((1+((SpecificHeatRatio-1)*pow(MachGhost,2)/2)),-1.0); 
	
	double Density = PressureGhost/(RealGasConstant*TemperatureGhost);

	double SoundSpeed = sqrt(SpecificHeatRatio*RealGasConstant*
		TemperatureGhost);

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

/*! \fn void BC(
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

\brief Function BC() implements the boundary condition. 
 Here ghost cell are used to implement the boundary condition. In simple 
words this function calculates the conserved variables for all ghost cells.
\param [in] ConservedVariables Pointer to the 4D vector where all 
the conserved variables of previous time step are stored.
\param [in] iFaceAreaVector Pointer to the 4D vector which has the 
area vector of all faces which are in "i" direction.  
\param [in] jFaceAreaVector Pointer to the 4D vector which has the 
area vector of all faces which are in "j" direction.  
\param [in] kFaceAreaVector Pointer to the 4D vector which has the 
area vector of all faces which are in "k" direction.  
\param [in,out] i0GhostConservedVariable Conserved variables in the ghost cells 
at i=0
\param [in,out] j0GhostConservedVariable Conserved variables in the ghost cells 
at j=0
\param [in,out] k0GhostConservedVariable Conserved variables in the ghost cells 
at k=0
\param [in,out] iNiGhostConservedVariable Conserved variables in the ghost cells 
at i=Ni
\param [in,out] jNjGhostConservedVariable Conserved variables in the ghost cells 
at j=Nj
\param [in,out] kNkGhostConservedVariable Conserved variables in the ghost cells 
at k=Nk
\param [in] Ni Number of cells in in "i" direction.  
\param [in] Nj Number of cells in in "j" direction.  
\param [in] Nk Number of cells in in "k" direction.  
\param [in] SpecificHeatRatio Specific heat ratio 
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
	/*! primitive variables at the inlet */
	double InletDensity ; 
	double InletXVelocity ; 
	double InletYVelocity ; 
	double InletZVelocity ; 
	double InletStaticPressure ;

	double ExitStaticPressure ;

	/*! Total quantities are given at the inlet */ 
	double InletTotalTemperature;	
	double InletTotalPressure;

	/*! \warning Inlet Mach will be given/used, 
	only in case of supersonic inlet*/
	double InletMach; 

	// Boundary condition options at all the boundaries 
	string BoundaryConditionati0 ;
	string BoundaryConditionatj0 ;
	string BoundaryConditionatk0 ;
	string BoundaryConditionatiNi ;
	string BoundaryConditionatjNj ;
	string BoundaryConditionatkNk ;

	// reading the input file 
	ifstream infile("inputfile");
	string aline;

	while(!infile.eof())// file ended
	{
		getline(infile,aline); // reading a line form file

		if (aline.find( "//" )!=0 && aline.empty()==false) 
		{
			if(aline.find("InletDensity")!=string::npos)
			{
				InletDensity = atof(aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletXVelocity")!=string::npos)
			{
				InletXVelocity = atof(aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletYVelocity")!=string::npos)
			{
				InletYVelocity = atof(aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletZVelocity")!=string::npos)
			{
				InletZVelocity = atof(aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletStaticPressure")!=string::npos)
			{
				InletStaticPressure = 
				atof (aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("ExitStaticPressure")!=string::npos)
			{
				ExitStaticPressure = 
				atof (aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletTotalTemperature")!=string::npos)
			{
				InletTotalTemperature = 
				atof (aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletTotalPressure")!=string::npos)
			{
				InletTotalPressure = 
				atof (aline.substr(aline.find("=")+1).c_str());
			}
			else if(aline.find("InletMach")!=string::npos)
			{
				InletMach = atof (aline.substr(aline.find("=")+1).c_str());
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
