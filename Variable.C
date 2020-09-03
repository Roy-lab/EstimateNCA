#include "Variable.H"

Variable::Variable()
{
}

Variable::~Variable()
{
}

int 
Variable::setName(const char* aStr)
{
	name.append(aStr);
	return 0;
}

int 
Variable::setName(const string& aStr)
{
	setName(aStr.c_str());
	return 0;
}

const string& 
Variable::getName()
{
	return name;
}

int
Variable::setID(int aId)
{
	vId=aId;
	return 0;
}

int
Variable::getID()
{
	return vId;
}
