#ifndef VARSEL_H
#define VARSEL_H

#include "Matrix.H"
#include <vector>
#include <algorithm>
using namespace std;

class VariableSelection
{
public:
	VariableSelection();
	~VariableSelection();
	int LASSO_MultiL(Matrix* X, Matrix* Y, vector<double>* ls, Matrix* betas);
	int LASSO_1L(Matrix* X, Matrix* Y, double l, Matrix* beta);
	int weightedLASSO_MultiL(Matrix* X, Matrix* Y, vector<double>* ls, Matrix* w, Matrix* betas);
	int weightedLASSO_1L(Matrix* X, Matrix* Y, double l, Matrix* w, Matrix* beta);
	int weightedLASSO_CV(Matrix* X, Matrix* Y, Matrix* w, Matrix* beta, int c);
private:
	double softTH(double v, double l);
	int testBetas(Matrix* X, Matrix* Y, Matrix* betas, vector<double>* errs);
	int getRandIndex(int* indices, int len);
};

#endif //VARSEL_H
