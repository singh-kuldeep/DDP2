/*! \file ArrayTester.h
    \brief Contains the function which checks the NaN/-NaN in array.	
*/
#ifndef ARRAYTESTER_H
#define ARRAYTESTER_H
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "iostream"
#include <vector>
#include "math.h"
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#include "Colortext.h"

using namespace std;

/*! \fn int test3DArray(string arrayname, vector<vector<vector<double> > > a, 
int Ni, int Nj, int Nk)
    \brief Checks the NaN in a 3D array.
    \param [in] arrayname Name of the array.
    \param [in] a Pointer to the array
    \param [in] Ni Array's first dimension 
    \param [in] Nj Array's second dimension
    \param [in] Nk Array's third dimension
    \return "0" if NaN and "1" if Not NaN 
*/
int test3DArray(string arrayname, vector<vector<vector<double> > > a, int Ni, 
	int Nj, int Nk)
{
	int flag=1;
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				if(isnan(a[i][j][k])==1)
				{
					cout << red(arrayname) << red("  has NaN at ") << i 
					<< ","<< j << "," << k << endl;
					return 0; 
				}
			}
		}
	}
	cout << blue(arrayname) << green("has 'NO' NaN parameter in it ") << endl;
	return 1;
}

/*! \fn int test4DArray(string arrayname, vector<vector<vector<vector<double> > > > a,
	int Ni, int Nj, int Nk, int Nl)
    \brief Checks the NaN in a 3D array.
    \param [in] arrayname Name of the array.
    \param [in] a Pointer to the array
    \param [in] Ni Array's first dimension 
    \param [in] Nj Array's second dimension
    \param [in] Nk Array's third dimension
    \param [in] Nl Array's fourth dimension
    \return "0" if NaN and "1" if Not NaN 
*/
int test4DArray(string arrayname, vector<vector<vector<vector<double> > > > a,
	int Ni, int Nj, int Nk, int Nl)
{
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				for (int l = 0; l < Nl; ++l)
				{
					if(isnan(a[i][j][k][l])==1)
					{
						cout << red(arrayname) << red(" has NaN at ") << i 
						<< ","<< j << "," << k << endl;
						return 0; 
					}
				}
			}
		}
	}
	
	cout << blue(arrayname) << green("has 'NO' NaN parameter in it ") << endl;
	return 1; 
}

/*! \fn int testConservedVariables(string arrayname, 
	vector<vector<vector<vector<double> > > > ConservedVariables, 
	int Ni, int Nj, int Nk, int Nl)
    \brief Checks whether density, pressure, energy are becoming NaN, 
    at any location.
    \param [in] arrayname Name of the array.
    \param [in] ConservedVariables Pointer to the array
    \param [in] Ni Array's first dimension 
    \param [in] Nj Array's second dimension
    \param [in] Nk Array's third dimension
    \param [in] Nl Array's fourth dimension
    \return "0" if NaN and "1" if Not NaN 
*/
int testConservedVariables(string arrayname, 
	vector<vector<vector<vector<double> > > > ConservedVariables, 
	int Ni, int Nj, int Nk, int Nl)
{
	double Density ; 
	/**\param Density Density at any location*/
	double XVelocity ;
	/**\param XVelocity X direction velocity at any location*/
	double YVelocity ;
	/**\param YVelocity Y direction velocity at any location*/
	double ZVelocity ;
	/**\param ZVelocity Z direction velocity at any location*/
	double pressure;
	/**\param pressure Pressure at any location*/
	double enenrgy ;
	/**\param energy Total energy at any location*/

	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				Density = ConservedVariables[i][j][k][0]; 
				XVelocity = ConservedVariables[i][j][k][1]/Density;
				YVelocity = ConservedVariables[i][j][k][2]/Density;
				ZVelocity = ConservedVariables[i][j][k][3]/Density;
				pressure = (1.4-1)*(ConservedVariables[i][j][k][4]-
				0.5*Density*(XVelocity*XVelocity + YVelocity*YVelocity +
				 ZVelocity*ZVelocity));
				enenrgy = ConservedVariables[i][j][k][4];

				if(isnan(Density)==1)
				{
					cout << red(arrayname) << red(" has NaN at ") 
					<< i << ","<< j << "," << k << red("because of Density") 
					<< endl;
					cout << " Density is " << Density << endl ;
					return 0; 
				}
				if(isnan(pressure)==1)
				{
					cout << red(arrayname) << red(" has NaN at ") 
					<< i << ","<< j << "," << k << red("because of pressure") 
					<< endl;
					cout << " pressure is " << pressure << endl ;
					return 0; 
				}
				if(isnan(enenrgy)==1)
				{
					cout << red(arrayname) << red(" has NaN at ") 
					<< i << ","<< j << "," << k << red("because of enenrgy") 
					<< endl;
					cout << " enenrgy is " << enenrgy << endl ;
					return 0; 
				}
			}
		}
	}
	
	cout << blue(arrayname) << green("has 'NO' NaN parameter in it ") << endl;
	return 1; 
}
#endif // array_tester.h ends here