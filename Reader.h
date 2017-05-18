/*! \file Reader.h
    \brief Contains the function ReadConservedVariables(), which reads the 
    conserved variables from the file.
  
    \date 18-May-2017 
*/

#ifndef READER_H
#define READER_H

#include <math.h>
#include "vector"
#include <iostream>
#include <fstream> 
#include <string> 
#include <cstdlib>
 
using namespace std;

/*! \fn void ReadConservedVariables(
	vector<vector<vector<vector<double> > > > & Uin,
	vector<vector<vector<vector<double> > > > & UinNew,
	int Ni, int Nj, int Nk)

\brief Reads the conserved variables from the file.
\param [in] Uin pointer to the conserved variable matrix of current time step
\param [in] UinNew pointer to the conserved variable matrix of the new time step
\param [in] Ni Number of cells in in "i" direction.  
\param [in] Nj Number of cells in in "j" direction.  
\param [in] Nk Number of cells in in "k" direction.
\return Magnitude of a 3D vector*/

void ReadConservedVariables(
	vector<vector<vector<vector<double> > > > & Uin,
	vector<vector<vector<vector<double> > > > & UinNew,
	int Ni, int Nj, int Nk)
{
	// reading the file
	ifstream ConservedQuantity("./Results/outputfiles/ConservedQuantity.csv");
   	for (int i = 0; i < Ni; ++i)
   	{
   		for (int j = 0; j < Nj; ++j)
   		{
		   string U;
		   getline(ConservedQuantity,U);
		   int length = U.size();
		   int l[4];
		   int commaPos = 0;

		   // noting the positions of ","
		   for (int Position = 0; Position < length; ++Position)
		   {
		   		if(U[Position]==',')
		   		{
		   			l[commaPos] = Position;
		   			commaPos++; 
		   		} 
		   }

   			for (int k = 0; k < Nk; ++k)
   			{
   				// converting string to float
				Uin[i][j][k][0] = atof(U.substr(0,(l[0]-0)).c_str());
				Uin[i][j][k][1] = atof(U.substr(l[0]+1,(l[1]-l[0]-1)).c_str());
				Uin[i][j][k][2] = atof(U.substr(l[1]+1,(l[2]-l[1]-1)).c_str());
				Uin[i][j][k][3] = atof(U.substr(l[2]+1,(l[3]-l[2]-1)).c_str());
				Uin[i][j][k][4] = atof(U.substr(l[3]+1,(l[4]-l[3]-1)).c_str());

				UinNew[i][j][k][0] = atof(U.substr(0,(l[0]-0)).c_str());
				UinNew[i][j][k][1] = 
				atof(U.substr(l[0]+1,(l[1]-l[0]-1)).c_str());
				UinNew[i][j][k][2] = 
				atof(U.substr(l[1]+1,(l[2]-l[1]-1)).c_str());
				UinNew[i][j][k][3] = 
				atof(U.substr(l[2]+1,(l[3]-l[2]-1)).c_str());
				UinNew[i][j][k][4] = 
				atof(U.substr(l[3]+1,(l[4]-l[3]-1)).c_str());
   			}
   		}
   	}
}
#endif // READER_H