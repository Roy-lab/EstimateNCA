#include "Graph.H"
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#include <gsl/gsl_rng.h>


Graph::Graph()
{
	regMngr = NULL;
	tgtMngr = NULL;
	adjMat  = NULL;
}

Graph::~Graph()
{
	if (adjMat != NULL)
	{
		delete adjMat;
	}
}

int 
Graph::setVariableManagers(VariableManager* rptr, VariableManager* tptr)
{
	regMngr = rptr;
	tgtMngr = tptr;
	map<int,Variable*> regVarSet = regMngr->getVariableSet();
	map<int,Variable*> tgtVarSet = tgtMngr->getVariableSet();
	int nRegs = regVarSet.size();
	int nTgts = tgtVarSet.size();
	//cout << nRegs <<"," << nTgts << endl;
	adjMat = new Matrix(nTgts,nRegs);
	return 0;
}

int
Graph::readNet(const char* inname)
{
	if (regMngr == NULL || tgtMngr == NULL)
	{
		cerr << "ERROR! (in readNet): VariableManagers are not set!" << endl;
		cerr << "Exiting..." << endl;
		exit(0);
	}

	ifstream inFile(inname);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string tfName;
		string tgName;
		double edgeStrength;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfName.append(tok);
			}
			else if(tokCnt==1)
			{
				tgName.append(tok);
			}
			else if(tokCnt==2)
			{
				edgeStrength=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		int uID = regMngr->getVarID(tfName.c_str());
		int vID = tgtMngr->getVarID(tgName.c_str());
		if (uID == -1 || vID == -1)
		{
			continue;
		}
		adjMat->setValue(edgeStrength,vID,uID);
	}
	inFile.close();
	return 0;
}

int
Graph::setNet(EdgeList* net)
{
	if (adjMat != NULL)
	{
		delete adjMat;
	}
	map<int,Variable*> regVarSet = regMngr->getVariableSet();
	map<int,Variable*> tgtVarSet = tgtMngr->getVariableSet();
	int nRegs = regVarSet.size();
	int nTgts = tgtVarSet.size();
	adjMat = new Matrix(nTgts,nRegs);
	NET_T* nnet = net->getNet();
	for (map<string,map<string,double>* >::iterator itr=nnet->begin(); itr!=nnet->end(); itr++)
	{
		string name = itr->first;
		map<string,double>* elist = itr->second;
		for (map<string,double>::iterator eitr=elist->begin(); eitr!=elist->end(); eitr++)
		{
			string gname = eitr->first;
			double v = eitr->second;
			int uID = regMngr->getVarID(name.c_str());
			int vID = tgtMngr->getVarID(gname.c_str());
			if (uID == -1 || vID == -1)
			{
				continue;
			}
			adjMat->setValue(v,vID,uID);
		}
	}
	return 0;
}

int
Graph::setAdj(Matrix* mat)
{
	if (adjMat != NULL)
	{
		delete adjMat;
	}
	adjMat = mat->copyMe();
	return 0;
}

Matrix*
Graph::getAdj()
{
	if (adjMat == NULL)
	{
		return NULL;
	}
	Matrix* temp = adjMat->copyMe();
	return temp;
}

int
Graph::removeRegs(vector<int>& ids)
{
	if (adjMat == NULL)
	{
		return 0;
	}
	adjMat->removeCols(ids);
	return 0;
}

int
Graph::removeTgts(vector<int>& ids)
{
	if (adjMat == NULL)
	{
		return 0;
	}
	adjMat->removeRows(ids);
	return 0;
}

int
Graph::randomizeWeights()
{
	if (adjMat == NULL)
	{
		return 0;
	}
	gsl_rng_env_setup();
	gsl_rng * r;
	r = gsl_rng_alloc (gsl_rng_default);
	gsl_rng_set(r, time(NULL));
	int nRegs = adjMat->getColCnt();
	int nTgts = adjMat->getRowCnt();
	for (int i=0;i<nTgts;i++)
	{
		for (int j=0;j<nRegs;j++)
		{
			double v = adjMat->getValue(i, j);
			if (v != 0)
			{
				v = gsl_rng_uniform (r);
				adjMat->setValue(v,i,j);
			}
		}
	}
	gsl_rng_free (r);
	return 0;
}

double
Graph::getEdgeWeight(int uID, int vID)
{
	if (adjMat == NULL)
	{
		return 0;
	}
	return adjMat->getValue(vID, uID);
}

int
Graph::writeNet(const char* oname)
{
	ofstream oFile(oname);
	int row = adjMat->getRowCnt();
	int col = adjMat->getColCnt();
	for (int i=0;i<row;i++)
	{
		Variable* vVar = tgtMngr->getVariableAt(i);
		for (int j=0;j<col;j++)
		{
			Variable* uVar = regMngr->getVariableAt(j);
			double v = adjMat->getValue(i,j);
			if (v != 0)
			{
				oFile << uVar->getName() << "\t" << vVar->getName() << "\t" << v << endl;
			}
		}
	}
	oFile.close();
	return 0;
}

EdgeList*
Graph::getNetMap()
{
	EdgeList* onet = new EdgeList;
	int row = adjMat->getRowCnt();
	int col = adjMat->getColCnt();
	for (int j=0;j<col;j++)
	{
		Variable* uVar = regMngr->getVariableAt(j);
		string tname = uVar->getName();
		map<string,double>* elist = new map<string,double>;
		for (int i=0;i<row;i++)
		{
			Variable* vVar = tgtMngr->getVariableAt(i);
			string gname = vVar->getName();
			double v = adjMat->getValue(i,j);
			if (v != 0)
			{
				onet->addEdge(tname,gname,v);
			}
		}
	}
	return onet;
}
