#ifndef _GRAPH
#define _GRAPH

#include "CommonTypes.H"
#include "VariableManager.H"
#include "Matrix.H"
#include "EdgeList.H"
#include <vector>
using namespace std;

class Graph
{
	public:
		Graph();
		~Graph();
		int setVariableManagers(VariableManager* rptr, VariableManager* tptr);
		int readNet(const char* inname);
		int setAdj(Matrix* mat);
		int setNet(EdgeList* net);
		Matrix* getAdj();
		int removeRegs(vector<int>& ids);
		int removeTgts(vector<int>& ids);
		int randomizeWeights();
		double getEdgeWeight(int uID, int vID);
		int writeNet(const char* oname);
		EdgeList* getNetMap();
	private:
		Matrix* adjMat;
		VariableManager* regMngr;
		VariableManager* tgtMngr;
};

#endif //_GRAPH
