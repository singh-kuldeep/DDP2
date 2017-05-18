/*! \file AllFacesFluxAUSM.h
    \brief Calculates the flux at all the cell faces in the domain using the 
    AUSM scheme. 

    \date 18-May-2017 
*/
#ifndef ALLFACESFLUXAUSM_h
#define ALLFACESFLUXAUSM_h

#include "vector"

#include "ArrayTester.h"
#include "NetFluxAUSM.h" //AUSM
#include "BC.h"

/*! \fn void flux(
	vector<vector<vector<vector<double> > > >	& iFacesFlux,
	vector<vector<vector<vector<double> > > >	& jFacesFlux,
	vector<vector<vector<vector<double> > > >	& kFacesFlux,
	vector<vector<vector<vector<double> > > >	iFaceAreaVector,
	vector<vector<vector<vector<double> > > >	jFaceAreaVector,
	vector<vector<vector<vector<double> > > >	kFaceAreaVector,
	vector<vector<vector<vector<double> > > >	ConservedVariables,
	int Ni, int Nj, int Nk,
	string gamma, double SpecificHeatRatio)
	
	\brief function flux() calculates the flux at every face in the domain
	\param [in,out] &iFacesFlux Pointer to an array which stores the flux at "i" 
	interfaces
	\param [in,out] &jFacesFlux Pointer to an array which stores the flux at "j" 
	interfaces
	\param [in,out] &kFacesFlux Pointer to an array which stores the flux at "k" 
	interfaces
	\param [in] &iFacesFlux An array which has area vectors of the "i" 
	interfaces
	\param [in] &jFacesFlux An array which has area vectors of the "j" 
	interfaces
	\param [in] &kFacesFlux An array which has area vectors of the "k" 
	interfaces
	*/
void flux(
	vector<vector<vector<vector<double> > > >	& iFacesFlux,
	vector<vector<vector<vector<double> > > >	& jFacesFlux,
	vector<vector<vector<vector<double> > > >	& kFacesFlux,
	vector<vector<vector<vector<double> > > >	iFaceAreaVector,
	vector<vector<vector<vector<double> > > >	jFaceAreaVector,
	vector<vector<vector<vector<double> > > >	kFaceAreaVector,
	vector<vector<vector<vector<double> > > >	ConservedVariables,
	int Ni, int Nj, int Nk,
	string gamma, double SpecificHeatRatio)
{
	// Creating a 4D vector object
	typedef vector<double> Dim1;
	typedef vector<Dim1> Dim2;
	typedef vector<Dim2> Dim3;
	typedef vector<Dim3> Dim4;

	std::vector<double> LeftConservedVariables(5);
	/**\param LeftConservedVariables Temporary vector which store the "left" 
	cell conserved variables*/ 
	std::vector<double> RightConservedVariables(5);
	/**\param RightConservedVariables Temporary vector which store the "right" 
	cell conserved variables*/

	/**\param i0GhostConservedVariable Ghost cell Conserved variables array 
	at i = 0*/
	/**\param j0GhostConservedVariable Ghost cell Conserved variables array 
	at j = 0*/
	/**\param k0GhostConservedVariable Ghost cell Conserved variables array 
	at k = 0*/
	/**\param iNiGhostConservedVariable Ghost cell Conserved variables array 
	at i = Ni*/
	/**\param jNjGhostConservedVariable Ghost cell Conserved variables array 
	at j = Nj*/
	/**\param kNkGhostConservedVariable Ghost cell Conserved variables array 
	at k = Nk*/
	Dim4 i0GhostConservedVariable(1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	Dim4 iNiGhostConservedVariable(1,Dim3(Nj,Dim2(Nk,Dim1(5))));
	Dim4 j0GhostConservedVariable(Ni,Dim3(1,Dim2(Nk,Dim1(5))));
	Dim4 jNjGhostConservedVariable(Ni,Dim3(1,Dim2(Nk,Dim1(5))));
	Dim4 k0GhostConservedVariable(Ni,Dim3(Nj,Dim2(1,Dim1(5))));
	Dim4 kNkGhostConservedVariable(Ni,Dim3(Nj,Dim2(1,Dim1(5))));

	/*! Before every time step we need to have proper value in all the ghost 
		cells So, function BC() takes care of this*/
	BC(ConservedVariables,
	iFaceAreaVector,jFaceAreaVector,kFaceAreaVector,
	i0GhostConservedVariable,j0GhostConservedVariable,
	k0GhostConservedVariable,iNiGhostConservedVariable,
	jNjGhostConservedVariable,kNkGhostConservedVariable,
	Ni, Nj, Nk,SpecificHeatRatio);
	
	/*Checking whether BC() function has assigned correct values to the ghost 
	cells*/
	
	/*only for debugging purpose. Do not uncomment while running the solver*/  
	#if 0 
	if(test4DArray("i0GhostConservedVariable",i0GhostConservedVariable,
	1,Nj,Nk,5)==0)
	{
		// return 0;
	}
	if(test4DArray("iNiGhostConservedVariable",iNiGhostConservedVariable,
	1,Nj,Nk,5)==0)
	{
		// return 0;
	}
	if(test4DArray("j0GhostConservedVariable",j0GhostConservedVariable,
	Ni,1,Nk,5)==0)
	{
		// return 0;
	} 
	if(test4DArray("jNjGhostConservedVariable",jNjGhostConservedVariable,
	Ni,1,Nk,5)==0)
	{
		// return 0;
	}
	if(test4DArray("k0GhostConservedVariable",k0GhostConservedVariable,
	Ni,Nj,1,5)==0)
	{
		// return 0;
	}
	if(test4DArray("kNkGhostConservedVariable",kNkGhostConservedVariable,
	Ni,Nj,1,5)==0)
	{
		// return 0;
	}
	#endif

	// i faces flux calculation
	for (int i = 0; i < Ni+1; ++i)
	{

		for (int j = 0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{	
				// i interface volume
				if(i == 0)
				{ 
					LeftConservedVariables = i0GhostConservedVariable[0][j][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the -0.5 interface 
					netfluxAUSM irightface(LeftConservedVariables,
					RightConservedVariables,
					iFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);


					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(irightface.NetFlux[l])==1)
						{
							cout << "irightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						iFacesFlux[i][j][k][l] = irightface.NetFlux[l]; 
					}
				}
				else if(i == Ni)
				{
					LeftConservedVariables = ConservedVariables[i-1][j][k];
					RightConservedVariables =iNiGhostConservedVariable[0][j][k];
					
					// Calculating the flux at the Ni-0.5 interface 
					netfluxAUSM irightface(LeftConservedVariables,
					RightConservedVariables,
					iFaceAreaVector[i][j][k], gamma,SpecificHeatRatio);
						
					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(irightface.NetFlux[l])==1)
						{
							cout << "irightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						iFacesFlux[i][j][k][l] = irightface.NetFlux[l]; 
					}
					
				}
				else
				{
					LeftConservedVariables = ConservedVariables[i-1][j][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the i-0.5 interface 
					netfluxAUSM irightface(LeftConservedVariables,
					RightConservedVariables,
					iFaceAreaVector[i][j][k], gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(irightface.NetFlux[l])==1)
						{
							cout << "irightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						iFacesFlux[i][j][k][l] = irightface.NetFlux[l]; 
					}

				}
			}
		}
	}

	// j faces flux calculation
	for (int i = 0; i < Ni; ++i)
	{
		for (int j =0; j < Nj+1; ++j)
		{
			for (int k = 0; k < Nk; ++k)
			{	
				if(j == 0)
				{
					LeftConservedVariables = j0GhostConservedVariable[i][0][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the -0.5 interface 
					netfluxAUSM jrightface(LeftConservedVariables,
					RightConservedVariables,
					jFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(jrightface.NetFlux[l])==1)
						{
							cout << "jrightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						jFacesFlux[i][j][k][l] = jrightface.NetFlux[l]; 
					}
				}
				else if(j == Nj)
				{
					LeftConservedVariables = ConservedVariables[i][j-1][k];
					RightConservedVariables = jNjGhostConservedVariable[i][0][k];

					// Calculating the flux at the Nj-0.5 interface 
					netfluxAUSM jrightface(LeftConservedVariables,
					RightConservedVariables,
					jFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(jrightface.NetFlux[l])==1)
						{
							cout << "jrightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						jFacesFlux[i][j][k][l] = jrightface.NetFlux[l]; 
					}
				}
				else
				{
					LeftConservedVariables = ConservedVariables[i][j-1][k];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the j-0.5 interface 
					netfluxAUSM jrightface(LeftConservedVariables,
					RightConservedVariables,
					jFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(jrightface.NetFlux[l])==1)
						{
							cout << "jrightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						jFacesFlux[i][j][k][l] = jrightface.NetFlux[l]; 
					}
				}

			}
		}
	}
	
	// k faces flux calculation
	for (int i = 0; i < Ni; ++i)
	{
		for (int j =0; j < Nj; ++j)
		{
			for (int k = 0; k < Nk+1; ++k)
			{
				if(k == 0)
				{						
					LeftConservedVariables = k0GhostConservedVariable[i][j][0];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the -0.5 interface 
					netfluxAUSM krightface(LeftConservedVariables,
					RightConservedVariables,
					kFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(jrightface.NetFlux[l])==1)
						{
							cout << "jrightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						kFacesFlux[i][j][k][l] = krightface.NetFlux[l]; 
					}
				}
				else if(k == Nk)
				{
					LeftConservedVariables = ConservedVariables[i][j][k-1];
					RightConservedVariables = kNkGhostConservedVariable[i][j][0];

					// Calculating the flux at the Nk-0.5 interface 
					netfluxAUSM krightface(LeftConservedVariables,
					RightConservedVariables,
					kFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(krightface.NetFlux[l])==1)
						{
							cout << "krightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						kFacesFlux[i][j][k][l] = krightface.NetFlux[l]; 
					}
				}
				else
				{
					LeftConservedVariables = ConservedVariables[i][j][k-1];
					RightConservedVariables = ConservedVariables[i][j][k];

					// Calculating the flux at the j-0.5 interface 
					netfluxAUSM krightface(LeftConservedVariables,
					RightConservedVariables,
					kFaceAreaVector[i][j][k],gamma,SpecificHeatRatio);

					for (int l = 0; l < 5; ++l)
					{
						#if 0 // Mainly for code testing 
						if(isnan(krightface.NetFlux[l])==1)
						{
							cout << "krightface.NetFlux["<<l<<"] is NaN at[" 
							<< i << "," << j << " " << k <<"]" << endl;
							cout << "Check the line number "<< __LINE__ <<
							" in TVDMainSolver.cpp" << endl;
							break;
						}
						#endif
						kFacesFlux[i][j][k][l] = krightface.NetFlux[l]; 
					}
				}
			}
		}
	}
}
#endif // ALLFACESFLUXAUSM_h ends here