#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "prng.h"
#include "montecarlo.h"
#include "hmc.h"
#include "MatrixMan.h"

/*
 * add \frac{\delta S[pmat]}{\delta pmat_{\sigma}} to pmom_{\sigma} where \sigma=1,2,3
 */
int deltaX3MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps)
{
	int i,j,a,b,c;
	doublecomplex *pcomm[nummat], *tmpmat1, *tmpmat2, *tmpmat[nummat];
	double tr=0;

	for(i=0;i<nummat;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmpmat[i], 0, sizeof(doublecomplex)*size*size);
	}
	tmpmat1 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
	tmpmat2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);

	Comm(pcomm[0], pmat[0], pmat[1], size, 1.0);				//Calculating commutators
	Comm(pcomm[1], pmat[1], pmat[2], size, 1.0);
	Comm(pcomm[2], pmat[2], pmat[0], size, 1.0);

	for(a=0;a<nummat;a++)
	{
		if(a==0)			//labels for different matrices in the array
		{
			b=1;			//value for second entry in commutator
			c=2;			//value for first entry in commutator
		}
		else if(a==1)
		{
			b=2;
			c=0;
		}
		else
		{
			b=0;
			c=1;
		}

		memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
		memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);

		Comm3(tmpmat1, pmat[b], pcomm[a], size, 1.0);
		Comm3(tmpmat2, pmat[c], pcomm[c], size, 1.0);

		AddMat(tmpmat[a], tmpmat1, pcomm[b], size, eps*size*1.0, -eps*coupling*size);
		AddMat(tmpmat[a], tmpmat2, pcomm[b], size, -eps*size*1.0, -eps*coupling*size);
	}

	// implementation of tracelessness
	for(j=0;j<nummat;j++)
	{
		tr=0;
		for(i=0;i<(size-1);i++)
		{
			(tmpmat[j]+i*size+i)->r -= (tmpmat[j]+(size-1)*size+(size-1))->r;
			tr += (tmpmat[j]+i*size+i)->r;
		}
		(tmpmat[j]+(size-1)*size+(size-1))->r = -tr;
	}
	for(i=0;i<nummat;i++)
	{
		AddMat6(pmom[i], tmpmat[i], 1.0, size);
	}

	for(i=0;i<nummat;i++)
	{
		free(pcomm[i]);
		free(tmpmat[i]);
	}
	free(tmpmat1);
	free(tmpmat2);

	return 0;
}
