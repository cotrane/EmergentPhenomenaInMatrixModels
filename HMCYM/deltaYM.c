#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"
#include "MatrixMan.h"

/*
 * compute \frac{ \delta S}{\delta X_{\sigma} and add it to momentum matrix pmom_{\sigma}
 */
int deltaYM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double eps, double mass)
{
	int i,j,k,a,b,pos1,pos2,numcomm=0;
	double con, tr=0;
	doublecomplex **pcomm, *tmpmat[nummat];

	//initialization
	for(i=0;i<nummat;i++)
	{
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmpmat[i], 0, sizeof(doublecomplex)*size*size);
	}

	for(i=0;i<nummat;i++)
		numcomm += (nummat-1-i);

	pcomm = (doublecomplex**) malloc(sizeof(doublecomplex*)*numcomm);
	for(i=0;i<numcomm;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	//computes commutator
	for(a=0; a<(nummat-1); a++)
	{
		pos1=0;
		for(i=0; i<a; i++)
		{
			pos1+=((nummat-1)-i);
		}
		pos2=0;
		for(b=(a+1); b<nummat; b++)
		{
			Comm(pcomm[pos1+pos2], pmat[a], pmat[b], size, 1.0);
			pos2++;
		}
	}

	if(mass==0)
	{
		// computes YM Term
		for(a=0;a<nummat;a++)
		{
			for(b=0;b<nummat;b++)
			{
				if(b!=a)
				{
					if(a<b)
					{
						con=1.0;
						pos1=0;
						for(i=0; i<a; i++)
						{
							pos1+=((nummat-1)-i);
						}
						Commend2(tmpmat[a], pmat[b], pcomm[pos1+b-1-a], size, eps*size*con);
					}
					if(b<a)
					{
						con=-1.0;
						pos1=0;
						for(i=0; i<b; i++)
						{
							pos1+=((nummat-1)-i);
						}
						Commend2(tmpmat[a], pmat[b], pcomm[pos1+a-1-b], size, eps*size*con);
					}

				}
			}
		}
	}

	//mass term

	if(mass!=0)
	{
		// computes YM Term
		for(a=0;a<nummat;a++)
		{
			for(b=0;b<nummat;b++)
			{
				if(b!=a)
				{
					if(a<b)
					{
						con=1.0;
						pos1=0;
						for(i=0; i<a; i++)
						{
							pos1+=((nummat-1)-i);
						}
						Commend2(tmpmat[a], pmat[b], pcomm[pos1+b-1-a], size, eps*size*con/(mass*mass));
					}
					if(b<a)
					{
						con=-1.0;
						pos1=0;
						for(i=0; i<b; i++)
						{
							pos1+=((nummat-1)-i);
						}
						Commend2(tmpmat[a], pmat[b], pcomm[pos1+a-1-b], size, eps*size*con/(mass*mass));
					}

				}
			}
		}

		// computes variation of mass term and adds it to tmpmat
		for(i=0;i<nummat;i++)
		{
			for(j=0;j<size;j++)
			{
				for(k=(j+1);k<size;k++)
				{
					(tmpmat[i]+(j*size + k))->r += -2.*eps*size* (pmat[i]+(k*size + j))->r;
					(tmpmat[i]+(j*size + k))->i += -2.*eps*size* (pmat[i]+(k*size + j))->i;

					(tmpmat[i]+(k*size + j))->r = (tmpmat[i]+(j*size + k))->r;
					(tmpmat[i]+(k*size + j))->i = -(tmpmat[i]+(j*size + k))->i;
				}
				(tmpmat[i]+(j*size + j))->r += -2.*eps*size* (pmat[i]+(j*size + j))->r;
			}
		}
	}

	// make matrices traceless by subtracting sum of first (N-1) diagonal elements from N'th diagonal element
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
	// adding result to pmom
	for(i=0;i<nummat;i++)
	{
		AddMat6(pmom[i], tmpmat[i], 1.0, size);
	}

	for(i=0;i<nummat;i++)
	{
		free(tmpmat[i]);
	}
	for(i=0;i<numcomm;i++)
	{
		free(pcomm[i]);
	}
	free(pcomm);

	return 0;
}
