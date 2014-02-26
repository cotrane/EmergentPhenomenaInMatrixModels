#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"

int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size)
{
	int i,j,k;

	for(i=0;i<nummat;i++)
	{
		for(j=0;j<size;j++)
		{
			for(k=(j+1);k<size;k++)
			{
				(pmat[i]+(j*size + k))->r += eps * ((pmom[i]+(k*size + j))->r);
				(pmat[i]+(j*size + k))->i += eps * (pmom[i]+(k*size + j))->i;

				(pmat[i]+(k*size + j))->r = (pmat[i]+(j*size + k))->r;
				(pmat[i]+(k*size + j))->i = -(pmat[i]+(j*size + k))->i;
			}
			(pmat[i]+(j*size + j))->r += eps * (pmom[i]+(j*size + j))->r;
		}
	}

	return 0;
}

int addmom(doublecomplex *pmom[], doublecomplex *pmat[], double eps, int nummat, int size, double mass)
{
	int j,k,l;

	for(j=0;j<nummat;j++)
	{
		for(k=0;k<size;k++)
		{
			for(l=(k+1);l<size;l++)
			{
				(pmom[j]+(k*size + l))->r += -2*size*eps * (pmat[j]+(l*size + k))->r;
				(pmom[j]+(k*size + l))->i += -2*size*eps * (pmat[j]+(l*size + k))->i;

				(pmom[j]+(l*size + k))->r = (pmom[j]+(k*size + l))->r;
				(pmom[j]+(l*size + k))->i = -(pmom[j]+(k*size + l))->i;
			}
			(pmom[j]+(k*size + k))->r += -2*size*eps * (pmat[j]+(k*size + k))->r;
		}
	}
	return 0;
}
double mom(doublecomplex *pmom[], int nummat, int size)
{
	int i;
	double mom=0;

	for(i=0;i<nummat;i++)
	{
		diagMulti(&mom, pmom[i], pmom[i], size, 0.5);
	}

	return mom;
}

double action(doublecomplex *pmat[], int nummat, int size, double mass)
{
	int i;
	double S=0;

	for(i=0;i<nummat;i++)
	{
		diagMulti(&S, pmat[i], pmat[i], size, size);
	}

	return S;
}
