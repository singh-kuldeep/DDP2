// this function will go through the whole array and search for the NaN/-NaN

#include "iostream"
#include <vector>
#include "math.h"
using namespace std;

// 3d array
int test3DArray(string arrayname, vector<vector<vector<double> > > a, int Ni, int Nj, int Nk)
{
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				if(isnan(a[i][j][k])==1)
				{
					cout << arrayname << "  has NaN at " << i << ","<< j << "," << k << endl;
					return 0; 
				}
			}
		}
	}
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
						cout << arrayname << " has NaN at " << i << ","<< j << "," << k << endl;
						return 0; 
					}
				}
			}
		}
	}
}