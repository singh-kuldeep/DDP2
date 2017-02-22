// checked
#include "math.h"
#include "iostream"
#include "interface.h"

#define gamma 1.4
#define gasconstant 287.14
#define heatcapacityconstantvolume 717.5

using namespace std ;
// [The link text](http://example.net/ "Link title")
/*! 
 *  \brief     Pretty nice class.
 *  \details   This class is used to demonstrate a number of section commands.
 *  \author    John Doe
 *  \author    Jan Doe
 *  \version   4.1a
 *  \date      1990-2011
 *  \pre       First initialize the system.
 *  \bug       Not all memory is freed when deleting an object of this class.
 *  \warning   Improper use can crash your application
 *  \copyright GNU Public License.
 */
class diffusionfluxinterface
{
	public:
	double diffusionfluxvector[5] ;

	diffusionfluxinterface(vector<double>& vectorleftminus, vector<double>& vectorleft, vector<double>& vectorright,
	vector<double>& vectorrightplus, vector<double>& areavectorleft, vector<double>& areavectorright, vector<double>& areavectorrightplus,
	double volumeleftmins, double volumeleft, double volumeright, double volumerightplus, double deltat)
	{
		double thetai[5]; /*!< Detailed description after the member */
		double thetaiplus[5]; /**< enum value 1 */
		double gvector[5] ;
		double gvectorplus[5] ;
		double signale ;
		double betainterface[5] ;
		double phiinterface[5] ;
		double Zvectorinterface[5];
		double pshiinterface[5] ;
	
	interface left(vectorleftminus,vectorleft,areavectorleft,volumeleftmins,volumeleft,deltat) ;
	interface right(vectorleft, vectorright, areavectorright,volumeleft,volumeright, deltat) ;
	interface rightplus(vectorright,vectorrightplus,areavectorrightplus,volumeright,volumerightplus,deltat);

//! A pure virtual member.
    /*!
      \sa testMe()
      \param c1 the first argument.
      \param c2 the second argument.
    */

/*! 
*  A list of events:
*    - mouse events
*         -#  The distance between \f$(x_1,y_1)\f$ and \f$(x_2,y_2)\f$ is \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2}\f$.
*         -# mouse move event
*         -# mouse click event\n
*            More info about the click event.
*         -# mouse double click event
*    - keyboard events
*         1. <a href="./www.google.com">link text</a> 
*         2. key up event
*         2. <a href="linkURL">link text</a> 
*
*
*  More text here.
*/
/*!	
* <table>
* <caption id="multi_row">Complex table</caption>
* <tr><th>Column 1                      <th>Column 2        <th>Column 3
* <tr><td rowspan="2">cell row=1+2,col=1<td>cell row=1,col=2<td>cell row=1,col=3
* <tr><td rowspan="2">cell row=2+3,col=2                    <td>cell row=2,col=3
* <tr><td>cell row=3,col=1                                  <td rowspan="2">cell row=3+4,col=3
* <tr><td colspan="2">cell row=4,col=1+2
* <tr><td>cell row=5,col=1              <td colspan="2">cell row=5,col=2+3
* <tr><td colspan="2" rowspan="2">cell row=6+7,col=1+2      <td>cell row=6,col=3
* <tr>                                                      <td>cell row=7,col=3
* <tr><td>cell row=8,col=1              <td>cell row=8,col=2\n
*   <table>
*     <tr><td>Inner cell row=1,col=1<td>Inner cell row=1,col=2
*     <tr><td>Inner cell row=2,col=1<td>Inner cell row=2,col=2
*   </table>
*   <td>cell row=8,col=3
*   <ul>
*     <li>Item 1
*     <li>Item 2
*   </ul>
* </table>	
*/

/*! 	
*/
// here I am somehow trying to defie the gi and gi+1
for (int i = 0; i < 5; ++i)
{
	if (right.gvectorinterface[i] >= 0)
	{
		signale = 1.0;
	}
	else
	{
		signale = -1.0 ;
	}
	// this is pura jugaad so change karna hai baad me "gvector[i]" ka expression
	double tempB = (left.gvectorinterface[i]*signale)  ;
	double tempA = fabs(right.gvectorinterface[i]) ;
	tempA = min(tempA,tempB) ;
	double temp_zero=0.0;
	tempB = max(temp_zero, tempA) ;
	gvector[i] = signale *  tempB ;

	// gvector[i] = signale * max(0.0, min( fabs(right.gvectorinterface[i]),
	// 		left.gvectorinterface[i]*signale ) ) ;
}

 for (int i = 0; i < 5; ++i)
{
	if (rightplus.gvectorinterface[i] >= 0)
	{
		signale = 1.0;
	}
	else
	{
		signale = -1.0 ;
	}
// this is pura jugaad so change karna hai baad me "gvector[i]" ka expression
	double tempA = fabs(rightplus.gvectorinterface[i]) ;
	double tempB = (right.gvectorinterface[i]*signale)  ;
	tempA = min(tempA,tempB) ;
	double temp_zero=0.0;
	tempB = max(temp_zero, tempA) ;
	gvectorplus[i] = signale * tempB ;
	// gvector[i] = signale * max(0.0, min( fabs(right.gvectorinterface[i]),
	// 		left.gvectorinterface[i]*signale ) ) ;
	
 }
// theta i term defying 
for (int i = 0; i < 5; ++i)
	{
		// here I have doubt about not equal to symbol because it cant be exect 0.00000 so most of the 
		//time we endup choosing theta i = 0.0
		if((fabs(right.alphavectorinterface[i]) + fabs(left.alphavectorinterface[i]))!=0.0)
		{
			thetai[i] = fabs(right.alphavectorinterface[i] - left.alphavectorinterface[i])/
				(fabs(right.alphavectorinterface[i]) + fabs(left.alphavectorinterface[i])) ;
		}
		else
		{
			thetai[i] = 0.0;
		}
	}

// theta i plus 
for (int i = 0; i < 5; ++i)
{
	// here I have doubt about not equal to symbol because it cant be exect 0.00000 so most of the
	//time we endup choosing theta i = 0.0
	if((fabs(rightplus.alphavectorinterface[i])+fabs(right.alphavectorinterface[i]))!=0.0)
	{
		thetaiplus[i]=fabs(rightplus.alphavectorinterface[i]-right.alphavectorinterface[i])/
			(fabs(rightplus.alphavectorinterface[i]) + fabs(right.alphavectorinterface[i])) ;
	}
	else
	{
		thetaiplus[i] = 0.0 ;
	}
}

// beta defying
// before that we nee to define the omega values which are constant 
double omega[5] = {0.25, 1.0, 1.0, 1.0, 0.25} ;
for (int i = 0; i < 5; ++i)
{
	betainterface[i] = 1.0 + omega[i]*(max(thetai[i],thetaiplus[i]));
}
/////////////////////////////////
// phi defying
for (int i = 0; i < 5; ++i)
{
	if (right.alphavectorinterface[i] != 0)
	{
		phiinterface[i] = (gvectorplus[i] - gvector[i]) / right.alphavectorinterface[i] ;
	}
	else
	{
		phiinterface[i] = 0.0 ;
	}
}

// now re-defying the Zinterface
for (int i = 0; i < 5; ++i)
{
	Zvectorinterface[i] = right.Zvectorinterface[i] + betainterface[i]*phiinterface[i] ;
}

// pshiinterface re-defying (for invisid wall flow)
for (int i = 0; i < 5; ++i)
{
	pshiinterface[i] = pow(Zvectorinterface[i],2) + 0.25 ;
}

// re-defying the pshiinterface
// // again we need to specify the deltaf = 0.2
// 	double deltaf = 0.2 ;
// 	for (int i = 0; i < 5; ++i)
// 	{
// 		if (Zvectorinterface[i] >= deltaf)
// 		{
// 			pshiinterface[i] = fabs(Zvectorinterface[i]) ;
// 		}
// 		else
// 		{
// 			pshiinterface[i] = 0.5*(pow(Zvectorinterface[i],2) + pow(deltaf,2))/ deltaf ;
// 		}
// 	}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// finally calculating the interface diffusion term = diffusionfluxvetor
// before that we can calculate tempvector as
double tempvector[5];
for (int i = 0; i < 5; ++i)
{
	tempvector[i] = ( betainterface[i]*(gvector[i]+gvectorplus[i]) -
	 pshiinterface[i]*right.alphavectorinterface[i] )/deltat ;
}
// now the diffusionfluxvector
for (int i = 0; i < 5; ++i)
{
	diffusionfluxvector[i] = 0.0 ;
	for (int l = 0; l < 5; ++l)
	{
		diffusionfluxvector[i] += right.eigenvectormatrix[i][l] * tempvector[l] ;
	}
}
};
	// ~diffusionfluxinterface();
	
};