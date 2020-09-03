#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include "Variable.H"
#include "Error.H"

#include <string>
#include <map>
#include <vector>

using namespace std;

class VariableManager
{
	public:
		VariableManager();
		~VariableManager();
		Error::ErrorCode readVariables(const char*);
		int getVarID(const char*);
		int getVarID(const string&);
		map<int,Variable*>& getVariableSet();
		Variable* getVariableAt(int);

		//Merge two variable managers
		int mergeVarSets(VariableManager* varMngr);
		//Remove variables
		int removeVarsByID(vector<int>& ids);

		VariableManager* copyMe();
		VariableManager* copyMeWithSuffix(const string& str);
		int getVarCnt();
	private:
		map<string,int> varNameIDMap;
		map<int,Variable*> variableSet;
};

#endif //_VARIABLE_MANAGER
