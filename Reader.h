
#include <iostream>
#include <fstream> 
#include <string> 
#include <vector>
using namespace std;

// int main()
void ReadConservedVariables(
	vector<vector<vector<vector<double> > > > & Uin,
	vector<vector<vector<vector<double> > > > & UinNew,
	int Ni, int Nj, int Nk)
{
	// std::vector<double> Uin(5);
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
		   for (int Position = 0; Position < length; ++Position)
		   {
		   		if(U[Position]==',')
		   		{
		   			l[commaPos] = Position;
		   			commaPos++; 
		   		} 
		   }
		   // double U0 = stof(U.substr(0,(l[0]-0)));
		   // double U1 = stof(U.substr(l[0]+1,(l[1]-l[0]-1)));
		   // double U2 = stof(U.substr(l[1]+1,(l[2]-l[1]-1)));
		   // double U3 = stof(U.substr(l[2]+1,(l[3]-l[2]-1)));
		   // double U4 = stof(U.substr(l[3]+1,(l[4]-l[3]-1)));
   			for (int k = 0; k < Nk; ++k)
   			{

			   Uin[i][j][k][0] = stof(U.substr(0,(l[0]-0)));
			   Uin[i][j][k][1] = stof(U.substr(l[0]+1,(l[1]-l[0]-1)));
			   Uin[i][j][k][2] = stof(U.substr(l[1]+1,(l[2]-l[1]-1)));
			   Uin[i][j][k][3] = stof(U.substr(l[2]+1,(l[3]-l[2]-1)));
			   Uin[i][j][k][4] = stof(U.substr(l[3]+1,(l[4]-l[3]-1)));

			   UinNew[i][j][k][0] = stof(U.substr(0,(l[0]-0)));
			   UinNew[i][j][k][1] = stof(U.substr(l[0]+1,(l[1]-l[0]-1)));
			   UinNew[i][j][k][2] = stof(U.substr(l[1]+1,(l[2]-l[1]-1)));
			   UinNew[i][j][k][3] = stof(U.substr(l[2]+1,(l[3]-l[2]-1)));
			   UinNew[i][j][k][4] = stof(U.substr(l[3]+1,(l[4]-l[3]-1)));
			   
			   // Uin[0] = stof(U.substr(0,(l[0]-0)));
			   // Uin[1] = stof(U.substr(l[0]+1,(l[1]-l[0]-1)));
			   // Uin[2] = stof(U.substr(l[1]+1,(l[2]-l[1]-1)));
			   // Uin[3] = stof(U.substr(l[2]+1,(l[3]-l[2]-1)));
			   // Uin[4] = stof(U.substr(l[3]+1,(l[4]-l[3]-1)));
			   

			   // cout << Uin[i][j][k][0]  << "  " << endl; //
			   // cout << i << ", " << j << ", " << k << endl; 
			   // cout << Uin[0] << " " << Uin[1] << " " <<
			   // Uin[2] << " " << Uin[3] << " " << Uin[4] << endl;
   			}
   		}
   	}
}