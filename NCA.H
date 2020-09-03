#ifndef _NCA
#define _NCA

#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include "Matrix.H"
#include "VariableManager.H"
#include "EvidenceManager.H"
#include "Graph.H"
#include "VariableSelection.H"

class NCA
{
	public:
		NCA(VariableManager* rptr, VariableManager* gptr, EvidenceManager* dptr, Graph* pptr, double l, int cv);
		~NCA();
		int estimate();

		int getData(VariableManager*& rptr, VariableManager*& gptr, Graph*& aptr, EvidenceManager*& tptr);
	private:
		VariableSelection* varsel;

		VariableManager* regMngr;
		VariableManager* tgtMngr;
		EvidenceManager* expMngr;
		EvidenceManager* tfaMngr;
		Graph* priorNet;
		Graph* regNet;
		int nGenes; //#Genes
		int nRegs; //#TFs
		int nSamples; //#Samples
		double lambda;
		int CV;

		map<int, vector<int>* > indeg;

		int doOnePrune();
		int findRed(Matrix* A1, vector<int>& rows, vector<int>& cols, int m1, int n1);
		int estimateP();
		int estimateA();
		int setMap();
};

#endif //_NCA
