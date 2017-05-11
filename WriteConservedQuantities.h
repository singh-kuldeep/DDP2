#ifndef WRITECONSERVEDQUANTITIES_H
#define WRITECONSERVEDQUANTITIES_H
void WriteConserveredQuantities(
vector<vector<vector<vector<double> > > > ConservedVariables,
int Ni, int Nj, int Nk)
{
	// storing the all conserved variables in one plane
	ofstream kullu_2D ;
	kullu_2D.open("./Results/outputfiles/ConservedQuantity.csv");
	#if 1
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
	#endif

	#if 0
	for (int i = 0; i < Ni+2; ++i)
	{
		for (int j = 0; j < Nj+2; ++j)
		{
			if(i==0 && j==0)
			{
				kullu_2D << ConservedVariables[i][j][Nk/2][0] << "," << 
				ConservedVariables[i][j][Nk/2][1] <<","<< 
				ConservedVariables[i][j][Nk/2][2] << "," <<
				ConservedVariables[i][j][Nk/2][3] << "," <<
				ConservedVariables[i][j][Nk/2][4] << endl ;
			}
			else if(i==0 && j==Nj+1)
			{
				kullu_2D << ConservedVariables[i][Nj-2][Nk/2][0] << "," << 
				ConservedVariables[i][Nj-2][Nk/2][1] <<","<< 
				ConservedVariables[i][Nj-2][Nk/2][2] << "," <<
				ConservedVariables[i][Nj-2][Nk/2][3] << "," <<
				ConservedVariables[i][Nj-2][Nk/2][4] << endl ;
			}
			else if(i==Ni+1 && j==0)
			{
				kullu_2D << ConservedVariables[Ni-2][j][Nk/2][0] << "," << 
				ConservedVariables[Ni-2][j][Nk/2][1] <<","<< 
				ConservedVariables[Ni-2][j][Nk/2][2] << "," <<
				ConservedVariables[Ni-2][j][Nk/2][3] << "," <<
				ConservedVariables[Ni-2][j][Nk/2][4] << endl ;
			}
			else if(i==Ni+1 && j==Nj+1)
			{
				kullu_2D << ConservedVariables[Ni-2][Nj-2][Nk/2][0] << "," << 
				ConservedVariables[Ni-2][Nj-2][Nk/2][1] <<","<< 
				ConservedVariables[Ni-2][Nj-2][Nk/2][2] << "," <<
				ConservedVariables[Ni-2][Nj-2][Nk/2][3] << "," <<
				ConservedVariables[Ni-2][Nj-2][Nk/2][4] << endl ;
			}
			else if(i==0 && j!=0 && j!=Nj+1)
			{
				kullu_2D << ConservedVariables[i][j-1][Nk/2][0] << "," << 
				ConservedVariables[i][j-1][Nk/2][1] <<","<< 
				ConservedVariables[i][j-1][Nk/2][2] << "," <<
				ConservedVariables[i][j-1][Nk/2][3] << "," <<
				ConservedVariables[i][j-1][Nk/2][4] << endl ;
			}
			else if(i==Ni+1 && j!=0 && j!=Nj+1)
			{
				kullu_2D << ConservedVariables[i-2][j-1][Nk/2][0] << "," << 
				ConservedVariables[i-2][j-1][Nk/2][1] <<","<< 
				ConservedVariables[i-2][j-1][Nk/2][2] << "," <<
				ConservedVariables[i-2][j-1][Nk/2][3] << "," <<
				ConservedVariables[i-2][j-1][Nk/2][4] << endl ;
			}
			else if(j==0 && i!=0 && i!=Ni+1)
			{
				kullu_2D << ConservedVariables[i-1][j][Nk/2][0] << "," << 
				ConservedVariables[i-1][j][Nk/2][1] <<","<< 
				ConservedVariables[i-1][j][Nk/2][2] << "," <<
				ConservedVariables[i-1][j][Nk/2][3] << "," <<
				ConservedVariables[i-1][j][Nk/2][4] << endl ;
			}
			else if(j==Nj+1 && i!=0 && i!=Ni+1)
			{
				kullu_2D << ConservedVariables[i-1][j-2][Nk/2][0] << "," << 
				ConservedVariables[i-1][j-2][Nk/2][1] <<","<< 
				ConservedVariables[i-1][j-2][Nk/2][2] << "," <<
				ConservedVariables[i-1][j-2][Nk/2][3] << "," <<
				ConservedVariables[i-1][j-2][Nk/2][4] << endl ;;
			}
			else if(i!=0 && i!=Ni+1 && j!=0 && j!=Nj+1)
			{
				kullu_2D << 
				ConservedVariables[i-1][j-1][Nk/2][0] <<","<< 
				ConservedVariables[i-1][j-1][Nk/2][1] <<","<< 
				ConservedVariables[i-1][j-1][Nk/2][2] <<","<<
				ConservedVariables[i-1][j-1][Nk/2][3] <<","<<
				ConservedVariables[i-1][j-1][Nk/2][4] <<endl ;
			}
		}
	}
	#endif   
}
#endif