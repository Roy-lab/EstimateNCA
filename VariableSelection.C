#include "VariableSelection.H"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

VariableSelection::VariableSelection()
{
}

VariableSelection::~VariableSelection()
{
}

int
VariableSelection::LASSO_MultiL(Matrix* X, Matrix* Y, vector<double>* ls, Matrix* betas)
{
	int m = X->getColCnt();
	sort(ls->begin(),ls->end(),greater<double>());
	betas->init(m,ls->size());
	Matrix* beta = new Matrix;
	beta->init(m,1);
	for (int i=0;i<ls->size();i++)
	{
		double l = (*ls)[i];
		//cout << "lambda: " << l << endl;
		LASSO_1L(X,Y,l,beta);
		for (int j=0;j<m;j++)
		{
			betas->setValue(beta->getValue(j,0),j,i);
		}
	}
	delete beta;
	return 0;
}

int
VariableSelection::LASSO_1L(Matrix* X, Matrix* Y, double l, Matrix* beta)
{
	double stopTH=0.01;
	int n = X->getColCnt();
	int m = X->getRowCnt();
	//beta->init(n,1);
	Matrix* oldBeta = new Matrix;
	oldBeta->init(n,1);
	//set old beta to inf
	oldBeta->setAllValues(10000);

	Matrix* residuals = new Matrix;
	Matrix* xb = new Matrix;
	residuals->init(m,1);
	xb->init(m,1);
	double e = 0;
	while (1)
	//for (int itr=0;itr<10;itr++)
	{
		for (int j=0;j<n;j++)
		{
			xb->setMultiplyMatrix(X,beta);
			//cout << "X*B" << endl;
			//xb->showMatrix();
			residuals->copyFrom(Y);
			residuals->subtractWithMatrix(xb);
			//cout << "Residuals" << endl;
			//residuals->showMatrix();
			Matrix* xj = X->getColMatrix(j);
			xj->dotMultiplyWithMatrix(residuals);
			double s = xj->sumCol(0);
			s = s/m;
			s = s+beta->getValue(j,0);
			//cout << "s: " << s << " ";
			s = softTH(s,l);
			//cout << ", soft: " << s << endl;
			beta->setValue(s,j,0);
			delete xj;
		}
		//e = (oldbeta - beta)' * (oldbeta - beta)
		oldBeta->subtractWithMatrix(beta);
		e = oldBeta->getFNorm();
		e = e/sqrt(n);
		//cout << "ERR: " << e << endl;
		if (e<stopTH)
			break;
		beta->copyTo(oldBeta);
	}
	delete oldBeta;
	delete residuals;
	delete xb;
	return 0;
}

int
VariableSelection::weightedLASSO_MultiL(Matrix* X, Matrix* Y, vector<double>* ls, Matrix* w, Matrix* betas)
{
	int m = X->getColCnt();
	sort(ls->begin(),ls->end(),greater<double>());
	betas->init(m,ls->size());
	Matrix* beta = new Matrix;
	beta->init(m,1);
	for (int i=0;i<ls->size();i++)
	{
		double l = (*ls)[i];
		//cout << "lambda: " << l << endl;
		weightedLASSO_1L(X,Y,l,w,beta);
		for (int j=0;j<m;j++)
		{
			betas->setValue(beta->getValue(j,0),j,i);
		}
	}
	delete beta;
	return 0;
}

int
VariableSelection::weightedLASSO_1L(Matrix* X, Matrix* Y, double l, Matrix* w, Matrix* beta)
{
	double stopTH=0.01;
	int n = X->getColCnt();
	int m = X->getRowCnt();
	//beta->init(n,1);
	Matrix* oldBeta = new Matrix;
	oldBeta->init(n,1);
	//set old beta to inf
	oldBeta->setAllValues(10000);

	Matrix* residuals = new Matrix;
	Matrix* xb = new Matrix;
	residuals->init(m,1);
	xb->init(m,1);
	double e = 0;
	while (1)
	//for (int itr=0;itr<10;itr++)
	{
		for (int j=0;j<n;j++)
		{
			double wj = w->getValue(j,0);
			xb->setMultiplyMatrix(X,beta);
			//cout << "X*B" << endl;
			//xb->showMatrix();
			residuals->copyFrom(Y);
			residuals->subtractWithMatrix(xb);
			//cout << "Residuals" << endl;
			//residuals->showMatrix();
			Matrix* xj = X->getColMatrix(j);
			xj->dotMultiplyWithMatrix(residuals);
			double s = xj->sumCol(0);
			s = s/m;
			s = s+beta->getValue(j,0);
			//cout << "s: " << s << " ";
			s = softTH(s,l*wj);
			//cout << ", soft: " << s << endl;
			beta->setValue(s,j,0);
			delete xj;
		}
		//e = (oldbeta - beta)' * (oldbeta - beta)
		oldBeta->subtractWithMatrix(beta);
		e = oldBeta->getFNorm();
		e = e/sqrt(n);
		//cout << "ERR: " << e << endl;
		if (e<stopTH)
			break;
		beta->copyTo(oldBeta);
	}
	delete oldBeta;
	delete residuals;
	delete xb;
	return 0;
}


double
VariableSelection::softTH(double v, double l)
{
	if (v>l)
		return v-l;
	if (v<-1*l)
		return v+l;
	return 0;
}

int
VariableSelection::weightedLASSO_CV(Matrix* X, Matrix* Y, Matrix* w, Matrix* beta, int c)
{
	int m = X->getRowCnt();
	//randomize the sample numbers
	int* indices = new int[m];
	getRandIndex(indices, m);
	//define borders of test sets
	int s = m/c;
	vector<int> borders;
	for (int i=0;i<=m;i+=s)
	{
		borders.push_back(i);
	}
	borders[borders.size()-1] = m;
	vector<double> ls;
	vector<double> allerrs;
	//define lambdas
	//for (int i=0;i<100;i++)
	for (int i=1;i<50;i++)
	{
		ls.push_back(double(i)/1000.0);
		allerrs.push_back(0);
	}
	sort(ls.begin(),ls.end(),greater<double>());

	//unsigned long int delta=0;
	//do cross validation
	for (int i=0;i<borders.size()-1;i++)
	{
		//unsigned long int t1;
		//unsigned long int t2;
		//struct rusage r_usage;
		//getrusage(RUSAGE_SELF,&r_usage);
		//t1 = r_usage.ru_maxrss;

		vector<int> trainids;
		vector<int> testids;
		for (int j=borders[i];j<borders[i+1];j++)
		{
			testids.push_back(indices[j]);
		}
		for (int j=borders[0];j<borders[i];j++)
		{
			trainids.push_back(indices[j]);
		}
		for (int j=borders[i+1];j<borders[borders.size()-1];j++)
		{
			trainids.push_back(indices[j]);
		}
		//remove test samples from training
		Matrix* Xtrain = X->copyMe();
		Xtrain->removeRows(testids);
		Matrix* Xtest  = X->copyMe();
		Xtest->removeRows(trainids);
		//remove training samples from test
		Matrix* Ytrain = Y->copyMe();
		Ytrain->removeRows(testids);
		Matrix* Ytest  = Y->copyMe();
		Ytest->removeRows(trainids);

		//test all betas
		Matrix* betas=new Matrix;
		weightedLASSO_MultiL(Xtrain,Ytrain,&ls,w,betas);
		//get the error on test
		vector<double>* errs = new vector<double>;
		testBetas(Xtest, Ytest, betas, errs);
		//cout << "errs->size(): " << errs->size() << endl;
		//sum the error over different folds of the same lambda
		for (int j=0;j<errs->size();j++)
		{
			double e0 = allerrs[j];
			double e  = errs->at(j);
			allerrs[j] = e+e0;
		}
		//cout << "BETAS!" << endl;
		//cout << "****************************" << endl;
		//betas->showMatrix();
		//cout << "****************************" << endl;

		delete Xtrain;
		delete Ytrain;
		delete Xtest;
		delete Ytest;
		delete betas;
		delete errs;
		
		//getrusage(RUSAGE_SELF,&r_usage);
		//t2 = r_usage.ru_maxrss;

		//delta = delta + (t2-t1);

		testids.clear();
		trainids.clear();
	}
	//printf("Delta = %ld\n",delta);
	//select best lambda with lowest error
	double l=ls[0];
	double e=allerrs[0];
	for (int j=0;j<allerrs.size();j++)
	{
		//cout << "j: " << j << ", e: " << allerrs[j] << ", lambda: " << ls[j] << endl;
		if (allerrs[j]<e)
		{
			e = allerrs[j];
			l = ls[j];
		}
	}
	//use the lambda on the whole data
	//cout << "selected l: " << l << endl;
	weightedLASSO_1L(X, Y, l, w, beta);
	delete indices;
	return 0;
}

int
VariableSelection::getRandIndex(int* indices, int len)
{
	for (int i=0;i<len;i++)
	{
		indices[i] = i;
	}
	gsl_rng_env_setup();
	gsl_rng * r;
	r = gsl_rng_alloc (gsl_rng_default);
	gsl_rng_set(r, time(NULL));
	gsl_ran_shuffle (r, indices, len, sizeof (int));
	gsl_rng_free (r);
	return 0;
}

int
VariableSelection::testBetas(Matrix* X, Matrix* Y, Matrix* betas, vector<double>* errs)
{
	errs->clear();
	Matrix* Pred = X->multiplyMatrix(betas);
	int n = Pred->getRowCnt();
	int l = Pred->getColCnt();
	//cout << "n: " << n << ", l: " << l << endl;
	double v1,v2;
	for (int i=0;i<l;i++)
	{
		double e=0;
		for (int j=0;j<n;j++)
		{
			v1 = Pred->getValue(j,i);
			v2 = Y->getValue(j,0);
			e = e + (v1-v2)*(v1-v2);
		}
		e = e/n;
		//cout << "e: " << e << endl;
		errs->push_back(e);
	}
	return 0;
}
