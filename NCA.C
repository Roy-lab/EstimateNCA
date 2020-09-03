#include "NCA.H"
#include <algorithm>

NCA::NCA(VariableManager* rptr, VariableManager* gptr, EvidenceManager* dptr, Graph* pptr, double l, int cv)
{
	varsel = new VariableSelection;
	lambda = l;
	CV = cv;

	//regMngr = rptr->copyMe();
	//tgtMngr = gptr->copyMe();
	regMngr = rptr;
	tgtMngr = gptr;

	//Matrix* E = dptr->getDataMat();
	//E->rowStandardize();
	//expMngr = new EvidenceManager;
	//expMngr->setVariableManager(tgtMngr);
	//expMngr->setDataMat(E);
	expMngr = dptr;

	nRegs = regMngr->getVarCnt();
	nGenes = tgtMngr->getVarCnt();

	Matrix* E = dptr->getDataMat();
	nSamples = E->getColCnt();

	Matrix* P = new Matrix(nRegs, nSamples);
	tfaMngr = new EvidenceManager();
	tfaMngr->setVariableManager(regMngr);
	tfaMngr->setDataMat(P);

	Matrix* W = pptr->getAdj();
	Matrix* A = W->copyMe();
	double maxVal = W->getMax();
	W->subtractScalar(maxVal);
	W->multiplyScalar(-1);
	priorNet = new Graph;
	priorNet->setVariableManagers(regMngr,tgtMngr);
	priorNet->setAdj(W);

	regNet = new Graph;
	regNet->setVariableManagers(regMngr,tgtMngr);
	regNet->setAdj(A);
	regNet->randomizeWeights();

	doOnePrune();

	delete E;
	delete P;
	delete W;
	delete A;
}

NCA::~NCA()
{
	delete varsel;
	//delete regMngr;
	//delete tgtMngr;
	//delete expMngr;
	delete tfaMngr;
	delete priorNet;
	delete regNet;
}

int
NCA::doOnePrune()
{
	Matrix* A = regNet->getAdj();
	vector<int> rows;
	vector<int> cols;
	findRed(A, rows, cols, nGenes, nRegs);
	delete A;

	regMngr->removeVarsByID(cols);
	tgtMngr->removeVarsByID(rows);
	tfaMngr->removeVarsByID(cols);
	expMngr->removeVarsByID(rows);
	priorNet->removeRegs(cols);
	priorNet->removeTgts(rows);
	regNet->removeRegs(cols);
	regNet->removeTgts(rows);

	nRegs = regMngr->getVarCnt();
	nGenes = tgtMngr->getVarCnt();

	setMap();

	return 0;
}

//Adjacency matrix A, assuming it is binary
//Rows to delete: rows
//Columns to delete: cols
//m rows (genes)
//n columns (TFs)
int
NCA::findRed(Matrix* A1, vector<int>& rows, vector<int>& cols, int m1, int n1)
{
	rows.clear();
	cols.clear();
	for (int i=n1-1;i>=0;i--)
	{
		//assuming binary matrix
		int d = A1->countColEq(i,0);
		if ((m1-d)<=2)
		{
			//fprintf(stdout,"TF %d has degree <= 2\n",i);
			cols.push_back(i);
			continue;
		}
		//For now, don't remove TFs that are subset of other TFs
		///* 
		for (int j=n1-1;j>=0;j--)
		{
			if (i==j)
				continue;
			if (find(cols.begin(),cols.end(),j) != cols.end()) // Has been removed before!
			{
				continue;
			}
			int flag=1;
			for (int k=0;k<m1;k++)
			{
				double vik = A1->getValue(k,i);
				if (vik==0)//TF i doesn't have an edge to gene k
					continue;
				double vjk = A1->getValue(k,j);
				if (vjk==0)//TF j doesn't have an edge to gene k
				{
					flag = 0;
					break;
				}
				//TF i and TF j both have an edge for gene k
				//If it happens for all genes, it means targets of i are subset of targets of j
				//we should remove i
			}
			if (flag)
			{
				//fprintf(stdout,"TF %d is subset of %d\n",i,j);
				cols.push_back(i);
				break;
			}
		}
		//*/
	}
	for (int i=0;i<m1;i++)
	{
		int flag = 1;
		for (int j=0;j<n1;j++)
		{
			double vij = A1->getValue(i,j);
			//TF j controls gene i, and TF j is not removed
			if (vij != 0 && find(cols.begin(),cols.end(),j)==cols.end())
			{
				flag = 0;
				break;
			}
		}
		if (flag)
		{
			rows.push_back(i);
			//fprintf(stdout,"Gene %d is only regulated by removed TFs\n",i);
		}
	}
	//double rc = gsl_blas_dasum(cols);
	//double rr = gsl_blas_dasum(rows);
	//fprintf(stdout,"Removed %d TFs and %d Genes\n",(int)rc,(int)rr);
	
	return 0;
}

int
NCA::estimateP()
{
	Matrix* A = regNet->getAdj();
	Matrix* E = expMngr->getDataMat();

	Matrix* AT   = A->transMatrix();
	Matrix* ATA  = AT->multiplyMatrix(A);
	Matrix* ATAI = ATA->invMatrix();
	Matrix* ATE  = AT->multiplyMatrix(E);
	Matrix* P    = ATAI->multiplyMatrix(ATE);
	P->rowStandardize();
	tfaMngr->setDataMat(P);
	
	delete P;
	delete A;
	delete E;
	delete AT;
	delete ATA;
	delete ATAI;
	delete ATE;
	return 0;
}

int
NCA::estimateA()
{
	Matrix* A = regNet->getAdj();
	Matrix* W = priorNet->getAdj();
	Matrix* E = expMngr->getDataMat();
	Matrix* P = tfaMngr->getDataMat();

	Matrix* Y = new Matrix(nSamples,1);

	for (int i=0;i<nGenes;i++)
	{

		vector<int>* Aindex = indeg[i];

		gsl_vector_view y = E->getRowView(i);
		Y->setCol(&y.vector,0);

		Matrix* P0 = new Matrix(nSamples,Aindex->size());
		Matrix* B  = new Matrix(Aindex->size(),1);
		Matrix* W0 = new Matrix(Aindex->size(),1);

		for (int j=0;j<Aindex->size();j++)
		{
			gsl_vector_view Pj = P->getRowView(Aindex->at(j));
			P0->setCol(&Pj.vector,j);
			W0->setValue(W->getValue(i,Aindex->at(j)),j,0);
		}

		//varsel->LASSO_1L(P0, Y, lambda, B);
		if (CV == 0)
		{
			varsel->weightedLASSO_1L(P0, Y, lambda, W0, B);
		}
		else
		{
			varsel->weightedLASSO_CV(P0, Y, W0, B, CV);
		}

		for (int j=0;j<Aindex->size();j++)
		{
			double v = B->getValue(j,0);
			A->setValue(v,i,Aindex->at(j));
		}
		delete B;
		delete P0;
		delete W0;
	}
	regNet->setAdj(A);
	delete A;
	delete W;
	delete E;
	delete P;
	return 0;
}

int
NCA::estimate()
{
	double v  = 1000;
	double v2 = 0;
	//for (int itr=0; itr<20; itr++)
	int itr=0;
	while (v>0.0001)
	{
		Matrix* oldP = tfaMngr->getDataMat();
		fprintf(stdout,"NCA, iter: %d, err: %.4f (%.4f)\n", itr,v,v2);
		cout << "#TFs:" << regMngr->getVarCnt() << endl;
		estimateP();
		estimateA();
		Matrix* P  = tfaMngr->getDataMat();
		Matrix* er = oldP->subtractMatrix(P);
		v2 = er->getFNorm();
		v = v2/(nRegs*nSamples);
		itr ++;

		doOnePrune();

		delete oldP;
		delete P;
		delete er;
	}
	fprintf(stdout,"NCA, iter: %d, err: %.4f (%.4f)\n", itr,v,v2);
	return 0;
}

int
NCA::setMap()
{
	double v=0;
	for (map<int, vector<int>* >::iterator itr=indeg.begin(); itr!=indeg.end(); itr++)
	{
		vector<int>* Aindex = itr->second;
		delete Aindex;
	}
	indeg.clear();
	for (int i=0;i<nGenes;i++)
	{
		vector<int>* Aindex = new vector<int>;
		for (int j=0;j<nRegs;j++)
		{
			v = regNet->getEdgeWeight(j,i); 
			if (v != 0)
				Aindex->push_back(j);
		}
		indeg[i] = Aindex;
	}
	return 0;
}

int
NCA::getData(VariableManager*& rptr, VariableManager*& gptr, Graph*& aptr, EvidenceManager*& tptr)
{
	rptr = regMngr->copyMe();
	gptr = tgtMngr->copyMe();

	aptr = new Graph;
	aptr->setVariableManagers(rptr,gptr);
	Matrix* A = regNet->getAdj();
	aptr->setAdj(A);
	delete A;

	tptr = new EvidenceManager;
	tptr->setVariableManager(rptr);
	Matrix* P = tfaMngr->getDataMat();
	tptr->setDataMat(P);
	delete P;

	return 0;
}
