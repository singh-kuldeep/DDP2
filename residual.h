#ifndef RESIDUAL_H
#define RESIDUAL_H
void residual(vector<double> & Residual,int iteration, vector<vector<vector<vector<double> > > > ConservedVariables,
vector<vector<vector<vector<double> > > > ConservedVariablesNew, int Ni, int Nj,
int Nk)
{
	// Residual calculation after each time step and writing the all 
	//residuals into the file
	/**\param DensityResidual Density residual*/
	/**\param xMomentumResidual x Momentum residual*/
	/**\param yMomentumResidual y Momentum residual*/
	/**\param zMomentumResidual z Momentum residual*/
	/**\param Energy residual */
	// double DensityResidual = 0.0 ; 
	// double xMomentumResidual = 0.0 ; 
	// double yMomentumResidual = 0.0 ; 
	// double zMomentumResidual = 0.0 ; 
	// double EnergyResidual = 0.0 ; 

	Residual[0] = 0;
	Residual[1] = 0;
	Residual[2] = 0;
	Residual[3] = 0;
	Residual[4] = 0;
	int TotalGridPoints = Ni*Nj*Nk ; 
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			// DensityResidual   += pow((ConservedVariablesNew[i][j][Nk/2][0] -
			// 	ConservedVariables[i][j][Nk/2][0]),2);
			// xMomentumResidual += pow((ConservedVariablesNew[i][j][Nk/2][1] - 
			// 	ConservedVariables[i][j][Nk/2][1]),2);     
			// yMomentumResidual += pow((ConservedVariablesNew[i][j][Nk/2][2] - 
			// 	ConservedVariables[i][j][Nk/2][2]),2);     
			// zMomentumResidual += pow((ConservedVariablesNew[i][j][Nk/2][3] - 
			// 	ConservedVariables[i][j][Nk/2][3]),2);     
			// EnergyResidual    += pow((ConservedVariablesNew[i][j][Nk/2][4] - 
			// 	ConservedVariables[i][j][Nk/2][4]),2);    
			
			Residual[0] += pow((ConservedVariablesNew[i][j][Nk/2][0] -
				ConservedVariables[i][j][Nk/2][0]),2);
			Residual[1] += pow((ConservedVariablesNew[i][j][Nk/2][1] - 
				ConservedVariables[i][j][Nk/2][1]),2);     
			Residual[2] += pow((ConservedVariablesNew[i][j][Nk/2][2] - 
				ConservedVariables[i][j][Nk/2][2]),2);     
			Residual[3] += pow((ConservedVariablesNew[i][j][Nk/2][3] - 
				ConservedVariables[i][j][Nk/2][3]),2);     
			Residual[4]    += pow((ConservedVariablesNew[i][j][Nk/2][4] - 
				ConservedVariables[i][j][Nk/2][4]),2);    
			


			// This is to stop the simulation automatically if nan occurs
			if(isnan(sqrt(Residual[0]))==1)
			{
				cout <<"sqrt(DensityResidual) is NaN at line "<<__LINE__<< endl;
			}
		}
	}
	// cout << "TotalGridPoints" << TotalGridPoints << endl ;
	// kullu_mass << iteration << ","<< sqrt(DensityResidual/
	// (Ni*Nj))  << "," << sqrt(xMomentumResidual/(Ni*Nj))
	// << "," << sqrt(yMomentumResidual/(Ni*Nj)) <<","<< sqrt(
	// zMomentumResidual/(Ni*Nj)) << "," << sqrt(EnergyResidual/
	// (Ni*Nj)) << endl ;			

	cout << sqrt(Residual[0]) << endl ;
}
#endif // residual.h ends here 
