// just to print the colored text
#include "iostream"

string fail()
{
	return "  \033[1;31m failed the test case \033[0m\n ";
}

string pass()
{
	return "  \033[1;32m passed the test case \033[0m\n ";
}

string blue(string inputstring)
{
	string output = "  \033[1;36m ";
	output.append(inputstring);
	output.append(" \033[0m");
	return output;
}

string red(string inputstring)
{
	string output = "  \033[1;31m ";
	output.append(inputstring);
	output.append(" \033[0m");
	return output;
}

string green(string inputstring)
{
	string output = "  \033[1;32m ";
	output.append(inputstring);
	output.append(" \033[0m");
	return output;
}