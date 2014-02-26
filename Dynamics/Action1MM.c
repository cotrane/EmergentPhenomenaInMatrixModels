#include<stdlib.h>
#include<f2c.h>
#include<string.h>
#include "montecarlo.h"
#include "MatrixMan.h"

/*
 * compute action for harmonic oscillator
 */
double action1MM(doublecomplex *pmat[], int nummat, int size)
{
	int i;
	double S=0;

	for(i=0;i<nummat;i++)
	{
		diagMulti(&S, pmat[i], pmat[i], size, 0.5);
	}

	return S;
}
