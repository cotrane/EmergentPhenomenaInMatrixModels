#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"

/*
 * adds momentum matrix to pmat for traceless matrices
 */
int addmat(double pmat[], double pmom[], double eps, int size)
{
	int i;
	double trace=0;

	for(i=0;i<(size-1);i++)
	{
		pmat[i] += eps * pmom[i];
		trace += pmat[i];
	}
	pmat[size-1] = -trace;

	return 0;
}

/*
 * adds momentum matrix pmom to pmat for matrices with nonzero trace
 */
int addmattrace(double pmat[], double pmom[], double eps, int size)
{
	int i;

	for(i=0;i<size;i++)
	{
		pmat[i] += eps * pmom[i];
	}

	return 0;
}

/*
 * computes 1/2 * p^2 for traceless matrices; last element is no independent degree of freedom here and thus it will
 * not be added to mom (traceless diagonal matrices only have N-1 degrees of freedom)
 */
double mom(double pmom[], int size)
{
	int i;
	double mom=0;

	for(i=0;i<(size-1);i++)
	{
		mom += 0.5*pmom[i]*pmom[i];
	}

	return mom;
}

/*
 * computes 1/2 * p^2 for matrices with nonzero trace
 */
double momtrace(double pmom[], int size)
{
	int i;
	double mom=0;

	for(i=0;i<(size);i++)
	{
		mom += 0.5*pmom[i]*pmom[i];
	}

	return mom;
}
