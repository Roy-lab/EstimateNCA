#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <algorithm>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"


VariableManager::VariableManager()
{
	variableSet.clear();
	varNameIDMap.clear();
}

VariableManager::~VariableManager()
{
	for (VSET_ITER itr=variableSet.begin();itr!=variableSet.end();itr++)
	{
		Variable* var=itr->second;
		delete var;
	}
	variableSet.clear();
	varNameIDMap.clear();
}

//Reads the schema of the variables

Error::ErrorCode
VariableManager::readVariables(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[400000];
	int nodeCnt=0;

	while(inFile.good())
	{
		inFile.getline(buffer,400000);

		if(strlen(buffer)<=0)
		{
			continue;
		}

		char* tok=strtok(buffer,"\t");
		string geneName(tok);

		Variable* var = new Variable;
		var->setID(nodeCnt);
		var->setName(tok);
		variableSet[nodeCnt]=var;
		
		varNameIDMap[geneName]=nodeCnt;

		nodeCnt++;
	}

	inFile.close();

	cout <<"Read information about " << variableSet.size() << " variables " << endl;

	return Error::SUCCESS;
}

int
VariableManager::getVarID(const string& varName)
{
	if(varNameIDMap.find(varName)==varNameIDMap.end())
	{
		return -1;
	}
	int vId=varNameIDMap[varName];
	return vId;
}

int
VariableManager::getVarID(const char* varName)
{
	string varKey(varName);
	if(varNameIDMap.find(varKey)==varNameIDMap.end())
	{
		return -1;
	}
	int vId=varNameIDMap[varKey];
	return vId;
}

map<int,Variable*>& 
VariableManager::getVariableSet()
{
	return variableSet;
}

Variable* 
VariableManager::getVariableAt(int vId)
{
	if(variableSet.find(vId)==variableSet.end())
	{
		cout << "Illegal variable id " << vId << endl;
		return NULL;
	}
	return variableSet[vId];
}

int 
VariableManager::mergeVarSets(VariableManager* varMngr)
{
	for (map<string,int>::iterator itr=varMngr->varNameIDMap.begin(); itr!=varMngr->varNameIDMap.end(); itr++)
	{
		string name = itr->first;
		//If it's not int the current map, add it
		if (varNameIDMap.find(name) == varNameIDMap.end())
		{
			int curSize = variableSet.size();
			Variable* var = new Variable;
			var->setID(curSize);
			var->setName(name.c_str());
			varNameIDMap[name] = curSize;
			variableSet[curSize] = var;
		}
	}
	return 0;
}

int 
VariableManager::removeVarsByID(vector<int>& ids)
{
	if (ids.size()==0)
		return 0;
	//If the IDs are not in the varSet, we are doing it wrong!
	//Exit
	for (int i=0;i<ids.size();i++)
	{
		if (variableSet.find(ids[i]) == variableSet.end())
		{
			cerr << "ERROR! (in removeVarsByID), id out of range! " << ids[i] << endl;
			cerr << "EXITING..." << endl;
			exit(0);
		}
	}
	
	//sort the IDs
	sort(ids.begin(),ids.end(),less<int>());
	
	//clear the names
	varNameIDMap.clear();
	
	map<int,Variable*> tempVarset;
	Variable* var;
	//set names from 0 to the first ID
	for (int j=0;j<ids[0];j++)
	{
		var = variableSet[j];
		varNameIDMap[var->getName()] = j;
		tempVarset[j] = var;
	}
	for (int i=0;i<ids.size();i++)
	{
		//delete ids[i]
		var = variableSet[ids[i]];
		delete var;
		
		int b=ids[i]+1;//beginning, one after current ID
		int e=0;       //end
		//If it is the last ID, set end to the last variable
		if (i==ids.size()-1)
		{
			e = variableSet.size();
		}
		else //otherwise set it to the next ID
		{
			e = ids[i+1];
		}
		//Copy variableSet (from b to e) to tempVarset and update varNameIDMap
		for (int j=b;j<e;j++)
		{
			var = variableSet[j];
			int curSize = tempVarset.size();
			var->setID(curSize);
			varNameIDMap[var->getName()] = curSize;
			tempVarset[curSize] = var;
		}
	}
	variableSet.clear();
	variableSet = tempVarset;
	return 0;
}

VariableManager*
VariableManager::copyMe()
{
	VariableManager* temp = new VariableManager;
	for (VSET_ITER itr=variableSet.begin(); itr!=variableSet.end(); itr++)
	{
		int vID = itr->first;
		Variable* mVar = itr->second;
		Variable* nVar = new Variable;
		nVar->setID(mVar->getID());
		nVar->setName(mVar->getName());
		temp->variableSet[vID] = nVar;
	}
	for (map<string,int>::iterator itr=varNameIDMap.begin(); itr!=varNameIDMap.end(); itr++)
	{
		temp->varNameIDMap[itr->first] = itr->second;
	}
	return temp;
}

VariableManager*
VariableManager::copyMeWithSuffix(const string& str)
{
	VariableManager* temp = new VariableManager;
	for (VSET_ITER itr=variableSet.begin(); itr!=variableSet.end(); itr++)
	{
		int vID = itr->first;
		Variable* mVar = itr->second;
		Variable* nVar = new Variable;
		nVar->setID(mVar->getID());
		string name = mVar->getName();
		name.append(str);
		nVar->setName(name);
		temp->variableSet[vID] = nVar;

		temp->varNameIDMap[name] = vID;
	}
	return temp;
}

int
VariableManager::getVarCnt()
{
	return variableSet.size();
}
