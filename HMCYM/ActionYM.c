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
 * compute value of YM Action with optional mass-term
 */
double actionYM(doublecomplex *pmat[], int size, int nummat, double mass)
{
	int i, numcomm=0,pos1,pos2,mat1,mat2;
	double trYM=0;
	doublecomplex **pcomm;

	// compute number of commutators
	for(i=0;i<nummat;i++)
		numcomm += (nummat-1-i);

	// initialization
	pcomm = (doublecomplex**) malloc(sizeof(doublecomplex*)*numcomm);
	for(i=0;i<numcomm;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	// compute commutators and save them in 'pcomm'
	for(mat1=0; mat1<(nummat-1); mat1++)
	{
		pos1=0;
		for(i=0; i<mat1; i++)
		{
			pos1+=((nummat-1)-i);
		}
		pos2=0;
		for(mat2=(mat1+1); mat2<nummat; mat2++)
		{
			Comm(pcomm[pos1+pos2], pmat[mat1], pmat[mat2], size, 1.0);
			pos2++;
		}
	}

	//YM-Term

	// square commutators and save sum in trYM
	if(mass==0)
	{
		for(i=0; i<numcomm; i++)															//Squaring them
		{
			diagMulti(&trYM, pcomm[i], pcomm[i], size, -size/2.0);
		}
	}

	//Mass-Term

	// if mass not zero: compute sum of mass terms and trace of commutator squared and add to trYM
	if(mass!=0)
	{
		for(i=0; i<numcomm; i++)															//Squaring them
		{
			diagMulti(&trYM, pcomm[i], pcomm[i], size, -size/(2.0*mass*mass));
		}

		for(i=0;i<nummat;i++)
		{
			diagMulti(&trYM, pmat[i], pmat[i], size, size);
		}
	}

	for(i=0;i<numcomm;i++)
	{
		free(pcomm[i]);
	}
	free(pcomm);

	return trYM;
}


