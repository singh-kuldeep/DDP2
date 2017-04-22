#include "iostream"
#include "fstream"
#include <string> /* For strings */
using namespace std;
int main()
{
	int TotalIteration;
	
	string Scheme ;// = "AUSM";// or "Roe"/
	string InletBC ;
	string InitialCondition ;
	string GeometryOption ;
	
	double deltat;
	double TemperatureFreestream;
	double PressureFreestream;
	double MachFreestream;

	ifstream infile("inputfile");
	string aline;


	while(!infile.eof())// file ended
	{
		getline(infile,aline);

		if (aline.find( "//" )!=0 && aline.empty()==false) 
		{
			size_t TotalIterationFound = aline.find("TotalIteration");
			
			// size_t SchemeFound = aline.find("Scheme");
			size_t InletBCFound = aline.find("InletBC");
			size_t InitialConditionFound = aline.find("InitialCondition");
			size_t GeometryOptionFoud = aline.find("GeometryOption");
			
			size_t deltatfound = aline.find("deltat");
			size_t TemperatureFreestreamFound = aline.find("TemperatureFreestream");
			size_t PressureFreestreamFound = aline.find("PressureFreestream");
			size_t MachFreestreamFound = aline.find("MachFreestream");
			
			if(TotalIterationFound!=string::npos)
			{
				TotalIteration = stoi (aline.substr(aline.find("=")+1));
			}
			else if (aline.find("Scheme")!=string::npos)
			{
				Scheme = aline.substr(aline.find("=")+2); 
			}
			else if (InletBCFound!=string::npos)
			{
				InletBC = aline.substr(aline.find("=")+2);
			}
			else if (InitialConditionFound!=string::npos)
			{
				InitialCondition = aline.substr(aline.find("=")+2);
			}
			else if (GeometryOptionFoud!=string::npos)
			{
				GeometryOption = aline.substr(aline.find("=")+2);
			}			
			else if(deltatfound!=string::npos)
			{
				deltat = stod (aline.substr(aline.find("=")+1));
			}
			else if(TemperatureFreestreamFound!=string::npos)
			{
				TemperatureFreestream = stod (aline.substr(aline.find("=")+1));
			}
			else if(PressureFreestreamFound!=string::npos)
			{
				PressureFreestream = stod (aline.substr(aline.find("=")+1));
			}
			else if(MachFreestreamFound!=string::npos)
			{
				MachFreestream = stod (aline.substr(aline.find("=")+1));
			}						
		}
	}
	// cout << MachFreestream << InitialCondition << deltat << endl;

	return 0;
} 