#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <map>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

#include <gsl/gsl_randist.h>

#include "Error.H"
#include "VariableManager.H"
#include "Matrix.H"
#include "CommonTypes.H"

using namespace std;

class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();

		//This associates the variable information with the variable id
		int setVariableManager(VariableManager*);
		
		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile_Continuous(const char*);

		int setDataMat(Matrix* mat);
		Matrix* getDataMat();

		int removeVarsByID(vector<int>& ids);

		int writeEvidence(const char* oname);
		int updateEvidence(EvidenceManager* eMngr);
	private:
		int readExpFromFile(const char* inname, map<string,vector<double>* >& expMap);

		VariableManager* vMgr;

		Matrix* dataMat;
};

#endif //_EVIDENCE_MANAGER
