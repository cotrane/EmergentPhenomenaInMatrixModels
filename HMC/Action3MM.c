#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"

/*
 * computes action for d=3 Yang-Mills-Myers model
 */

double action3MM(doublecomplex *pmat[], double *trYM1, double *trCS1, double coupling, int size, int nummat)
{
	int i;
	double trCS=0, trYM=0, trrand =0;
	doublecomplex *pcomm[nummat];

	for(i=0;i<nummat;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	//Calculating commutators and saving them in pcomm
	Comm(pcomm[0], pmat[0], pmat[1], size, 1);
	Comm(pcomm[1], pmat[1], pmat[2], size, 1);
	Comm(pcomm[2], pmat[2], pmat[0], size, 1);

	//YM-Term
	// squaring commutators and adding result to trYM
	diagMulti(&trYM, pcomm[0], pcomm[0], size, 1);
	diagMulti(&trYM, pcomm[1], pcomm[1], size, 1);
	diagMulti(&trYM, pcomm[2], pcomm[2], size, 1);

	//Myers-Term
	// computing Myers term and adding it to trCS
	diagMulti(&trCS, pcomm[0], pmat[2], size, 1);

	// multiplying with constants
	*trYM1 = size*(-1/2.0 *(trYM));
	*trCS1 = -size*coupling*(2.0 *(trCS));

	// adding terms to final result in trrand
	trrand = *trYM1 + *trCS1;

	for(i=0;i<nummat;i++)
	{
		free(pcomm[i]);
	}
	return trrand;
}
