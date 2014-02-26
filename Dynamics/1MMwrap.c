
#include<stdlib.h>
#include<f2c.h>
#include"hmc.h"


double action1MMwrap(doublecomplex *pmat[], double coupling, int size, int nummat, double C2, double C3, double LAMBDA, double *S0, double *S1)
{
	double S=0;

	S = action1MM(pmat, nummat, size);

	return S;
}

int deltaX1MMwrap(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps, double C2, double C3, double kappa)
{
	addmom1MM(pmat, pmom, eps, nummat, size);

	return 0;
}

