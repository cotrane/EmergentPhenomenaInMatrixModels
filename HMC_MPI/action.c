#include<stdlib.h>
#include<f2c.h>
#include<string.h>
#include "montecarlo.h"
#include "MatrixMan.h"


double action(doublecomplex *pmat[], int nummat, int size, double kappa)
{
	int i;
	double S=0;

	for(i=0;i<nummat;i++)
	{
		diagMulti(&S, pmat[i], pmat[i], size, 0.5);
	}

	return S;
}
