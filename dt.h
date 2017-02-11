// this function will calculate the delta_t for each local cell at every iteration
#include <vector>
#include <math.h>
#include <iostream>
using namespace std;
double dt(int i,int j, int k, vector<vector<vector<double> > > delta_s, vector<vector<vector<vector<double> > > > variablesvector){
	double gamma = 1.4;
	double CFL = 0.01;

	double q = sqrt(pow((variablesvector[i][j][k][1]/variablesvector[i][j][k][0]),2) +
	pow((variablesvector[i][j][k][2]/variablesvector[i][j][k][0]),2) +
	pow((variablesvector[i][j][k][3]/variablesvector[i][j][k][0]),2) );

	double p = (gamma-1)*(variablesvector[i][j][k][4] - 0.5*variablesvector[i][j][k][0]*q*q);
	double a = sqrt(gamma*p/variablesvector[i][j][k][0]);
	double delta_t = (CFL*delta_s[i][j][k])/(a+q);
	cout <<  "dt   i " << i << " j " << j << " k " << k << " dt = " << delta_t << endl;
	return delta_t;
} 