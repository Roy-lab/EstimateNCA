#ifndef _NCALEARNER
#define _NCALEARNER

#include "VariableManager.H"
#include "EvidenceManager.H"
#include "Graph.H"

#include <vector>
#include <map>

using namespace std;

class NCALearner
{
	public:
		NCALearner(VariableManager* rptr, VariableManager* gptr, EvidenceManager* dptr, Graph* pptr, double l, int cv, int nr, const char* oname);
		~NCALearner();
		int start();
		int getData(VariableManager*& rptr, EvidenceManager*& tptr, const char* suff);
	private:
		double lambda;
		int CV;
		int numRep;

		VariableManager* regMngr;
		VariableManager* tgtMngr;
		EvidenceManager* expMngr;
		Graph* priorNet;
		char outdir[1025];

		int countNames(vector<VariableManager*>* allVars, map<string,int>& counts);
		int getMissingVarIDs(VariableManager* vMngr, map<string,int>& counts, vector<int>& ids);
		int filterAllMissing(vector<VariableManager*>* allTFs, vector<VariableManager*>* allGenes, 
						vector<EvidenceManager*>* allP, vector<Graph*>* allA);
		int collapseReps(VariableManager*& oregMngr, VariableManager*& otgtMngr, EvidenceManager*& otfaMngr, Graph*& onet,
						vector<VariableManager*>* allTFs,vector<VariableManager*>* allGenes,vector<EvidenceManager*>* allP,
						vector<Graph*>* allA);
		int getStdCorr(Matrix* A, Matrix* B, Matrix*& C);

		VariableManager* outregMngr;
		VariableManager* outtgtMngr;
		EvidenceManager* outtfaMngr;
		Graph* outnet;
};

#endif //_NCALEARNER
