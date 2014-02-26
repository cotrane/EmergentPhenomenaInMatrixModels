#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"

/*
 * compute action for 3d YM-Myers model
 */
double action3MM(doublecomplex *pmat[], double coupling, int size, int nummat)
{
	int i;
	double trCS=0, trYM=0, trrand =0;
	doublecomplex *pcomm[nummat];

	for(i=0;i<nummat;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	//Initial Energy

	Comm(pcomm[0], pmat[0], pmat[1], size, 1);				//Calculating commutators
	Comm(pcomm[1], pmat[1], pmat[2], size, 1);
	Comm(pcomm[2], pmat[2], pmat[0], size, 1);

	//YM-Term

	diagMulti(&trYM, pcomm[0], pcomm[0], size, 1);			//Squaring them
	diagMulti(&trYM, pcomm[1], pcomm[1], size, 1);
	diagMulti(&trYM, pcomm[2], pcomm[2], size, 1);

	//CS-Term

	diagMulti(&trCS, pcomm[0], pmat[2], size, 1);

	//Initial Energy

	trYM = size*(-1/2.0 *(trYM));
	trCS = -size*coupling*(2.0 *(trCS));

	trrand = trYM + trCS;

	for(i=0;i<nummat;i++)
	{
		free(pcomm[i]);
	}
	return trrand;
}

