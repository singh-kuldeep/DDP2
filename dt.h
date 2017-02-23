/*! \file dt.h
    \brief This header file conditions the function TimeStep() which calculate 
    the local time step for each cell at every iteration.
    \author Kuldeep Singh
    \date 2016
    \sa grid()
    \bug Currently not using this, because  grid() is not calculating ds value. 
    So recheck this function as well after fixing the grid() function.
    \param [in] i,j,k Cell location for which TimeStep is to be calculated 
    \param [in] delta_s ds value of the cell for which TimeStep is to be 
    calculated
    \param [IN] ConservedVariables Conserved variables vector 
    \param CFL Courant–Friedrichs–Lewy number
    \param Pressure Static Pressure
    \param VelocityMagnitude Magnitude of the velocity
    \param VelocitySound Speed of sound
    \param [out] TimeStep Time step (dt)
	\return double
*/

#include <vector>
#include <math.h>
// #include <iostream>
using namespace std;
double TimeStep(
	int i,int j, int k,
	vector<vector<vector<double> > > delta_s,
	vector<vector<vector<vector<double> > > > ConservedVariables)
{
	double SpecificHeatRatio = 1.4;
	double CFL = 0.01;

	double VelocityMagnitude = sqrt(pow((ConservedVariables[i][j][k][1]/
	ConservedVariables[i][j][k][0]),2) +
	pow((ConservedVariables[i][j][k][2]/ConservedVariables[i][j][k][0]),2) +
	pow((ConservedVariables[i][j][k][3]/ConservedVariables[i][j][k][0]),2) );

	double Pressure = (SpecificHeatRatio-1)*(ConservedVariables[i][j][k][4] -
	0.5*ConservedVariables[i][j][k][0]*VelocityMagnitude*VelocityMagnitude);
	double VelocitySound = sqrt(SpecificHeatRatio*Pressure/
	ConservedVariables[i][j][k][0]);
	double TimeStep = (CFL*delta_s[i][j][k])/(VelocitySound+VelocityMagnitude);
	// cout <<  "TimeStep   i " << i << " j " << j << " k " << k << 
	// " TimeStep = " << TimeStep << endl;
	return TimeStep;
} 