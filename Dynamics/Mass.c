#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"
#include "MatrixMan.h"


double Mass(doublecomplex *pmat[], int nummat, int size, double constant)
{
	int i;
	double out;

	for(i=0;i<nummat;i++)
	{
		diagMulti(&out, pmat[i], pmat[i], size, constant);
	}

	return out;
}


int deltaMass(doublecomplex *pout[], doublecomplex *pmat[], int nummat, int size, double constant)
{
	int j,k,l;

	for(j=0;j<nummat;j++)
	{
		for(k=0;k<size;k++)
		{
			for(l=(k+1);l<size;l++)
			{
				(pout[j]+(k*size + l))->r += -constant * (pmat[j]+(l*size + k))->r;
				(pout[j]+(k*size + l))->i += -constant * (pmat[j]+(l*size + k))->i;

				(pout[j]+(l*size + k))->r = (pout[j]+(k*size + l))->r;
				(pout[j]+(l*size + k))->i = -(pout[j]+(k*size + l))->i;
			}
			(pout[j]+(k*size + k))->r += -constant * (pmat[j]+(k*size + k))->r;
		}
	}

	return 0;
}
