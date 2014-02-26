#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include "montecarlo.h"
#include "hmc.h"
#include<string.h>
#include "model.h"
#include "MatrixMan.h"

int time_rev(doublecomplex *pmom[], doublecomplex *pmat[], int size, double eps, double ksi, int steps, int nummat, double alphatilde, double C2, double C3, double lambda)
{
	int i,k,l;
	doublecomplex *pmatold[nummat], *pmomold[nummat];

	for(i=0;i<nummat;i++)
	{
		pmatold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmomold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);

		memset(pmomold[i], 0, sizeof(doublecomplex)*size*size);
		memset(pmatold[i], 0, sizeof(doublecomplex)*size*size);
	}
	for(i=0;i<nummat;i++)
	{
		memcpy(pmatold[i], pmat[i], sizeof(doublecomplex)*size*size);
		memcpy(pmomold[i], pmom[i], sizeof(doublecomplex)*size*size);
	}


	for(i=0;i<nummat;i++)
	{
		for(k=0;k<size;k++)
		{
			for(l=0;l<size;l++)
			{
				(pmomold[i]+(k*size + l))->r = -1.0 * (pmomold[i]+(k*size + l))->r;
				(pmomold[i]+(k*size + l))->i = -1.0 * (pmomold[i]+(k*size + l))->i;
			}
		}
	}

	for(i=0;i<steps;i++)
	{
		addmat(pmatold, pmomold, ksi*eps, nummat, size, lambda);
		deltaaction_function(pmatold, pmomold, nummat, size, alphatilde, eps/2.0, C2, C3, lambda);
		addmat(pmatold, pmomold, eps*(1-2*ksi), nummat, size, lambda);
		deltaaction_function(pmatold, pmomold, nummat, size, alphatilde, eps/2.0, C2, C3, lambda);
		addmat(pmatold, pmomold, ksi*eps, nummat, size, lambda);
	}

	for(i=0;i<nummat;i++)
	{
		for(k=0;k<size;k++)
		{
			for(l=0;l<size;l++)
			{
				(pmomold[i]+(k*size + l))->r = -1.0 * (pmomold[i]+(k*size + l))->r;
				(pmomold[i]+(k*size + l))->i = -1.0 * (pmomold[i]+(k*size + l))->i;
			}
		}
	}

	for(i=0;i<steps;i++)
	{
		addmat(pmatold, pmomold, ksi*eps, nummat, size, lambda);
		deltaaction_function(pmatold, pmomold, nummat, size, alphatilde, eps/2.0, C2, C3, lambda);
		addmat(pmatold, pmomold, eps*(1-2*ksi), nummat, size, lambda);
		deltaaction_function(pmatold, pmomold, nummat, size, alphatilde, eps/2.0, C2, C3, lambda);
		addmat(pmatold, pmomold, ksi*eps, nummat, size, lambda);
	}

	printf("after timereversal:\n");
	for(i=0;i<nummat;i++)
	{
		printf("pmom[%d]:\n", i);
//		printmat(pmom[i], size);
		printmatdiff(pmomold[i], pmom[i], size);
		printf("pmat[%d]:\n", i);
//		printmat(pmat[i], size);
		printmatdiff(pmatold[i], pmat[i], size);
	}

	for(i=0;i<nummat;i++)
	{
		free(pmatold[i]);
		free(pmomold[i]);
	}
	return 0;
}
