/*
 * action3MMwrap.c
 *
 *  Created on: 3 Aug 2011
 *      Author: tkaltenbrunner
 */

#include<stdlib.h>
#include<f2c.h>
#include"hmc.h"


double action3MMwrap(doublecomplex *pmat[], double coupling, int size, int nummat, double C2, double C3, double LAMBDA, double *S0, double *S1)
{
	double S=0;

	S = action3MM(pmat, coupling, size, nummat);

	return S;
}

int deltaX3MMwrap(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps, double C2, double C3, double LAMBDA)
{
	deltaX3MM(pmat, pmom, nummat, size, coupling, eps);

	return 0;
}
