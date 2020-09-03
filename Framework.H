#ifndef _FRAMEWORK
#define _FRAMEWORK

#include <vector>
#include <map>
#include <string>

#include "Matrix.H"
#include "NCALearner.H"
#include "VariableManager.H"
#include "EvidenceManager.H"
#include "Graph.H"

#include <map>
#include <vector>
using namespace std;

class Framework
{
public:
	Framework();
	~Framework();
	int init(int argc, char** argv);
	int start();
private:
	VariableManager* regMngr;
	VariableManager* tgtMngr;
	VariableManager* allGenesMngr;
	EvidenceManager* evidMngr;
	Graph* priorNet;
	EdgeList* pnet;
	double lambda; //LASSO penalty
	int CV;//CV
	int numNcaRep;
	char outdir[1025];

	int printHelp(char* name);

	int inferNCA(EdgeList* initNet);
};

#endif //_FRAMEWORK
