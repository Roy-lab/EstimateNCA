#ifndef _VARIABLE
#define _VARIABLE

#include <map>
#include <string>

using namespace std;

class Variable;

typedef map<int,Variable*> VSET;
typedef map<int,Variable*>::iterator VSET_ITER;

class Variable
{
	public:
		Variable();
		~Variable();
		//This is how the variable will be used throughout
		//the system
		int setID(int);
		int getID();
		//This is the name assigned to the variable by the user
		int setName(const char*);
		int setName(const string&);
		const string& getName();
	private:
		string name;
		int vId;
};

#endif //_VARIABLE
