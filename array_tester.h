#ifndef ARRAYTESTER_H
#define ARRAYTESTER_H
// this function will go through the whole array and search for the NaN/-NaN

#include "iostream"
#include <vector>
#include "math.h"
#include "colortext.h"
using namespace std;

// 3d array
int test3DArray(string arrayname, vector<vector<vector<double> > > a, int Ni, int Nj, int Nk)
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
					cout << red(arrayname) << red("  has NaN at ") << i << ","<< j << "," << k << endl;
					return 0; 
				}
			}
		}
	}
	cout << blue(arrayname) << green("has 'NO' NaN parameter in it ") << endl;
	return 1;
}

// 4d array
int test4DArray(string arrayname, vector<vector<vector<vector<double> > > > a, int Ni, int Nj, int Nk, int Nl)
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
						cout << red(arrayname) << red(" has NaN at ") << i << ","<< j << "," << k << endl;
						return 0; 
					}
				}
			}
		}
	}
	
	cout << blue(arrayname) << green("has 'NO' NaN parameter in it ") << endl;
	return 1; 
}

// Conserverd variables array checking
int testConservedVariables(string arrayname, 
	vector<vector<vector<vector<double> > > > ConservedVariables, 
	int Ni, int Nj, int Nk, int Nl)
{
	double Density ; // density in live cell
	double XVelocity ;
	double YVelocity ;
	double ZVelocity ;
	double pressure;
	double enenrgy ;
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				Density = ConservedVariables[i][j][k][0]; // density in live cell
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
					<< i << ","<< j << "," << k << red("because of Density") << endl;
					cout << " Density is " << Density << endl ;
					return 0; 
				}
				if(isnan(pressure)==1)
				{
					cout << red(arrayname) << red(" has NaN at ") 
					<< i << ","<< j << "," << k << red("because of pressure") << endl;
					cout << " pressure is " << pressure << endl ;
					return 0; 
				}
				if(isnan(enenrgy)==1)
				{
					cout << red(arrayname) << red(" has NaN at ") 
					<< i << ","<< j << "," << k << red("because of enenrgy") << endl;
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