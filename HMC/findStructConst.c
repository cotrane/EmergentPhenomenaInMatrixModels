/*
 * findStructConst.c
 *
 *  Created on: 2 Nov 2012
 *      Author: tkaltenbrunner
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<f2c.h>
#include "MatrixMan.h"
#include "montecarlo.h"
#include "RepsSUn.h"
#include "hmc.h"

/*
 * computes structure constant for lie=2 (3MM) or lie=3 (8MM)
 */
int findStructConst(structconst2 **StructConst, int lie)
{
	int i,j,k;
	int numbasis=0, numNonZero=0;
	double trace;
	doublecomplex **basis, *comm;
	structconst2 *tmpSC;


	numbasis = lie*lie-1;
	basis = (doublecomplex**) malloc(sizeof(doublecomplex*)*numbasis);
	for(i=0;i<numbasis;i++)
	{
		basis[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*lie*lie);
		memset(basis[i], 0, sizeof(doublecomplex)*lie*lie);
	}
	comm = (doublecomplex*) malloc(sizeof(doublecomplex)*lie*lie);

	// obtain lie generators of current lie group
	RepsSUn(basis, lie, lie, numbasis, 1.0);

	// compute structure constant and saving values in struct StructConst
	for(i=0;i<numbasis;i++)
	{
		for(j=(i+1);j<numbasis;j++)
		{
			memset(comm, 0, sizeof(doublecomplex)*lie*lie);
			Comm(comm, basis[i], basis[j], lie, 1.0);
			for(k=(j+1);k<numbasis;k++)
			{
				trace=0;
				diagMulti(&trace, comm, basis[k], lie, 2.0);
				if(trace)
				{
					if(numNonZero==(*StructConst)[0].numEntries)
					{
						tmpSC = realloc(*StructConst, sizeof(structconst2)*((*StructConst)[0].numEntries+1));
						*StructConst = tmpSC;
						(*StructConst)[0].numEntries++;
					}
					(*StructConst)[numNonZero].a = i;
					(*StructConst)[numNonZero].b = j;
					(*StructConst)[numNonZero].c = k;
					(*StructConst)[numNonZero].val = trace;
					numNonZero++;
				}
			}
		}
	}

	return 0;
}

