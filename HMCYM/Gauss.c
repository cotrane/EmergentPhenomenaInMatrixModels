#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"

/*
 * compute variation of momentum, p^2 / 2, which is added during integration step to pmat; for traceless matrices
 */
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size)
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
				(pmat[i]+(j*size + k))->r += eps * ((pmom[i]+(k*size + j))->r);
				(pmat[i]+(j*size + k))->i += eps * (pmom[i]+(k*size + j))->i;

				(pmat[i]+(k*size + j))->r = (pmat[i]+(j*size + k))->r;
				(pmat[i]+(k*size + j))->i = -(pmat[i]+(j*size + k))->i;
			}
			(pmat[i]+(j*size + j))->r += eps * (pmom[i]+(j*size + j))->r;
			trace += (pmat[i]+(j*size + j))->r;
		}
		(pmat[i]+(size-1)*size + (size-1))->r = -trace;
	}

	return 0;
}
/*
 * compute variation of momentum, p^2 / 2, which is added during integration step to pmat
 */
int addmattrace(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size)
{
	int i,j,k;
	double trace=0;

	for(i=0;i<nummat;i++)
	{
		trace=0;
		for(j=0;j<(size);j++)
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
/*
 * compute value of momentum for hermitian traceless matrix
 */
double mom(doublecomplex *pmom[], int nummat, int size)
{
	int i,j,k;
	double mom=0;

	for(i=0;i<nummat;i++)
	{
		for(j=0;j<(size-1);j++)
		{
			for(k=0;k<size;k++)
			{
				mom += 0.5*((pmom[i]+j*size+k)->r * (pmom[i]+k*size+j)->r - (pmom[i]+j*size+k)->i * (pmom[i]+k*size+j)->i);
			}
		}
		for(k=0;k<(size-1);k++)
		{
			mom += 0.5*((pmom[i]+(size-1)*size+k)->r * (pmom[i]+k*size+(size-1))->r - (pmom[i]+(size-1)*size+k)->i * (pmom[i]+k*size+(size-1))->i);
		}
	}

	return mom;
}
/*
 * compute value of momentum for hermitian matrix
 */
double momtrace(doublecomplex *pmom[], int nummat, int size)
{
	int i,j,k;
	double mom=0;

	for(i=0;i<nummat;i++)
	{
		for(j=0;j<(size);j++)
		{
			for(k=0;k<size;k++)
			{
				mom += 0.5*((pmom[i]+j*size+k)->r * (pmom[i]+k*size+j)->r - (pmom[i]+j*size+k)->i * (pmom[i]+k*size+j)->i);
			}
		}
	}

	return mom;
}
