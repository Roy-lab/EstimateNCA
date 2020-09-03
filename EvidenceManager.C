#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "EvidenceManager.H"


EvidenceManager::EvidenceManager()
{
	vMgr = NULL;
	dataMat = NULL;
}

EvidenceManager::~EvidenceManager()
{
	if (dataMat!=NULL)
	{
		delete dataMat;
	}
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

int
EvidenceManager::removeVarsByID(vector<int>& ids)
{
	if (dataMat==NULL)
	{
		return 0;
	}
	dataMat->removeRows(ids);
	return 0;
}

int
EvidenceManager::readExpFromFile(const char* inname, map<string,vector<double>* >& expMap)
{
	ifstream inFile(inname);
	char buffer[500001];
	while(inFile.good())
	{
		inFile.getline(buffer,500000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string geneName;
		vector<double>* vals = new vector<double>;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else 
			{
				//vals->push_back(stod(string(tok)));
				vals->push_back(atoi(tok));
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if (expMap.find(geneName)!=expMap.end())
		{
			cerr << "ERROR! (in readExpFromFile) Shouldn't happen, gene \"" << geneName << "\" was repeated." << endl;
		}
		expMap[geneName] = vals;
	}
	return 0;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Continuous(const char* inname)
{
	if (vMgr == NULL)
	{
		cerr << "ERROR! (in loadEvidenceFromFile_Continuous): VariableManager is not set!" << endl;
		cerr << "Exiting..." << endl;
		exit(0);
	}
	map<string,vector<double>* > expMap;
	readExpFromFile(inname, expMap);
	map<int,Variable*> vset = vMgr->getVariableSet();
	vector<int> remIDs;
	for (VSET_ITER itr=vset.begin(); itr!=vset.end(); itr++)
	{
		Variable* var = itr->second;
		if (expMap.find(var->getName()) == expMap.end())
		{
			remIDs.push_back(var->getID());
		}
	}
	vMgr->removeVarsByID(remIDs);
	vset = vMgr->getVariableSet();
	map<string,vector<double>* >::iterator itr = expMap.begin();
	int nGenes = vset.size();
	int nSampels = itr->second->size();
	dataMat = new Matrix(nGenes,nSampels);
	for (itr=expMap.begin(); itr!=expMap.end(); itr++)
	{
		int vID = vMgr->getVarID(itr->first.c_str());
		vector<double>* vals = itr->second;
		if (vID == -1)
		{
			delete vals;
			continue;
		}
		for (int j=0;j<nSampels;j++)
		{
			dataMat->setValue(vals->at(j),vID,j);
		}
		delete vals;
	}
	expMap.clear();
	return Error::SUCCESS;
}

int
EvidenceManager::setDataMat(Matrix* mat)
{
	if (dataMat != NULL)
	{
		delete dataMat;
	}
	dataMat = mat->copyMe();
	return 0;
}

Matrix*
EvidenceManager::getDataMat()
{
	if (dataMat == NULL)
	{
		return NULL;
	}
	Matrix* temp = dataMat->copyMe();
	return temp;
}

int
EvidenceManager::writeEvidence(const char* oname)
{
	ofstream oFile(oname);
	int row = dataMat->getRowCnt();
	int col = dataMat->getColCnt();
	for (int i=0;i<row;i++)
	{
		Variable* var = vMgr->getVariableAt(i);
		oFile << var->getName();
		for (int j=0;j<col;j++)
		{
			oFile << "\t" << dataMat->getValue(i,j);
		}
		oFile << endl;
	}
	oFile.close();
	return 0;
}

//It assumes that VariableManager is already updated
int 
EvidenceManager::updateEvidence(EvidenceManager* eMngr)
{
	int nGenes = dataMat->getRowCnt();
	int nSampels = dataMat->getColCnt();
	//Resize the dataMat if the size of VariableManager and dataMat don't match
	if (nGenes < vMgr->getVarCnt())
	{
		int oldNGenes = nGenes;
		nGenes = vMgr->getVarCnt();
		Matrix* temp = new Matrix(nGenes,nSampels);
		for (int i=0;i<oldNGenes;i++)
		{
			for (int j=0;j<nSampels;j++)
			{
				temp->setValue(dataMat->getValue(i,j),i,j);
			}
		}
		delete dataMat;
		dataMat = temp;
	}
	VSET their_vset = eMngr->vMgr->getVariableSet();
	for (VSET_ITER itr=their_vset.begin(); itr!=their_vset.end(); itr++)
	{
		Variable* their_var = itr->second;
		int their_i = their_var->getID();
		string name = their_var->getName();
		int our_i   = vMgr->getVarID(name);
		if (our_i == -1)
		{
			cerr << "ERROR! (updateEvidence), two VariableManagers don't match!" << endl;
			cerr << "Exiting..." << endl;
			exit(0);
		}
		for (int j=0;j<nSampels;j++)
		{
			dataMat->setValue(eMngr->dataMat->getValue(their_i,j),our_i,j);
		}
	}
	return 0;
}
