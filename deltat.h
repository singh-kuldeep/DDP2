#ifndef DELTA_H
#define DELTA_H
// this function calculates the time step
double getLocalDeltaT(vector<double> ConservedVariables, double MinimumDistance, double CFL)
{
	double deltat = 1000; // time step

	double Density ;
	double Pressure ;
	double Velocity ;
	double VelocitySound;

	Density = ConservedVariables[0];

	Pressure = (SpecificHeatRatio -1)*
	(ConservedVariables[4] - 0.5*
	(pow(ConservedVariables[1],2)+
	pow(ConservedVariables[2],2)+
	pow(ConservedVariables[3],2))/Density ) ;  

	Velocity = sqrt(pow(ConservedVariables[1],2)+
	pow(ConservedVariables[2],2)+
	pow(ConservedVariables[3],2))/Density;
	
	VelocitySound = sqrt(SpecificHeatRatio*Pressure/
	Density);
	if(deltat > (CFL*MinimumDistance)/
		(Velocity+VelocitySound))
	{
		deltat = (CFL*MinimumDistance)/
		(Velocity+VelocitySound);
	}
	return deltat;
}	


double getGlobalDeltaT(vector<vector<vector<vector<double> > > > ConservedVariables,
vector<vector<vector<double> > > MinimumDistance, double CFL, int Ni, int Nj, int Nk)
{
	double deltat = 1000.0;
	for (int i = 0; i < Ni; ++i)
	{
		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{
				double LocalDeltaT = getLocalDeltaT(ConservedVariables[i][j][k], MinimumDistance[i][j][k],CFL);
				if(deltat > LocalDeltaT)
				{
					deltat = LocalDeltaT;
				}
			}
		}
	}
}	
#endif // de;tat.h ends here