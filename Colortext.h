/*! \file Colortext.h
    \brief Contains the function which prints the colorful texts. 
    
    \date 18-May-2017 
*/
#ifndef COLCOR_TEXT 
#define COLCOR_TEXT

#include "iostream"
#include <string> 

using namespace std;
/*!\fn string fail()
\brief Print the warning in red color if some test case fails
\return string
*/
string fail()
{
	return "  \033[1;31m failed the test case \033[0m\n ";
}

/*!\fn string pass()
\brief Print the success in green color if some test case gets pass
\return string
*/
string pass()
{
	return "  \033[1;32m passed the test case \033[0m\n ";
}

/*!\fn string blue(string inputstring)
\brief Return the string after appending the syntax of blue color
\param [IN] inputstring Input string 
\return string
*/
string blue(string inputstring)
{
	string output = " \033[1;36m ";
	output.append(inputstring);
	output.append(" \033[0m");
	return output;
}

/*!\fn string red(string inputstring)
\brief Return the string after appending the syntax of red color
\param [IN] inputstring Input string 
\return string
*/
string red(string inputstring)
{
	string output = "  \033[1;31m ";
	output.append(inputstring);
	output.append(" \033[0m");
	return output;
}

/*!\fn string green(string inputstring)
\brief Return the string after appending the syntax of green color
\param [IN] inputstring Input string 
\return string
*/
string green(string inputstring)
{
	string output = "  \033[1;32m ";
	output.append(inputstring);
	output.append(" \033[0m");
	return output;
}
#endif /* !FILE_FOO_SEEN */