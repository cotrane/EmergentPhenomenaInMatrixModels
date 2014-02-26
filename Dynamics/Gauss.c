#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"

/*
 * add \frac{\delta 0.5*pmom^2}{\delta pmom} to pmat
 * mass changes the frequency
 */
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, double mass)
{
	int i,j,k;
	double trace=0;

	for(i=0;i<nummat;i++)
	{
		trace=0;
		for(j=0;j<(size-1);j++)
		{
			for(k=(j+1);k<size;k++)
			{
				(pmat[i]+(j*size + k))->r += eps/mass * ((pmom[i]+(k*size + j))->r);
				(pmat[i]+(j*size + k))->i += eps/mass * (pmom[i]+(k*size + j))->i;

				(pmat[i]+(k*size + j))->r = (pmat[i]+(j*size + k))->r;
				(pmat[i]+(k*size + j))->i = -(pmat[i]+(j*size + k))->i;
			}
			(pmat[i]+(j*size + j))->r += eps/mass * (pmom[i]+(j*size + j))->r;
			trace += (pmat[i]+(j*size + j))->r;
		}
		(pmat[i]+(size-1)*size + (size-1))->r = -trace;
	}

	return 0;
}

/*
 * add \frac{\delta pmom^2/(2*mass)}{\delta pmom} to pmat
 */
int addmat1MM(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, double mass)
{
	int i,j,k;

	for(i=0;i<nummat;i++)
	{
		for(j=0;j<size;j++)
		{
			for(k=(j+1);k<size;k++)
			{
				(pmat[i]+(j*size + k))->r += eps/mass * ((pmom[i]+(k*size + j))->r);
				(pmat[i]+(j*size + k))->i += eps/mass * (pmom[i]+(k*size + j))->i;

				(pmat[i]+(k*size + j))->r = (pmat[i]+(j*size + k))->r;
				(pmat[i]+(k*size + j))->i = -(pmat[i]+(j*size + k))->i;
			}
			(pmat[i]+(j*size + j))->r += eps/mass * (pmom[i]+(j*size + j))->r;
		}
	}

	return 0;
}
/*
 * add \frac{\delta S[pmat]}{\delta pmat} to pmom
 */
int addmom1MM(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size)
{
	int j,k,l;

	for(j=0;j<nummat;j++)
	{
		for(k=0;k<size;k++)
		{
			for(l=(k+1);l<size;l++)
			{
				(pmom[j]+(k*size + l))->r += -eps * (pmat[j]+(l*size + k))->r;
				(pmom[j]+(k*size + l))->i += -eps * (pmat[j]+(l*size + k))->i;

				(pmom[j]+(l*size + k))->r = (pmom[j]+(k*size + l))->r;
				(pmom[j]+(l*size + k))->i = -(pmom[j]+(k*size + l))->i;
			}
			(pmom[j]+(k*size + k))->r += -eps * (pmat[j]+(k*size + k))->r;
		}
	}
	return 0;
}

/*
 * compute Tr(0.5 * pmom^2 / mass)
 */
double mom(doublecomplex *pmom[], int nummat, int size, double mass)
{
	int i;
	double mom=0;

	for(i=0;i<nummat;i++)
	{
		diagMulti(&mom, pmom[i], pmom[i], size, 0.5/mass);
	}

	return mom;
}
