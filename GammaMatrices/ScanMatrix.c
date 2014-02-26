#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>

/*
 * helper functions for generating \gamma-matrices
 * only needed for test of correctness
 */

int ScanMatrixZero(doublecomplex *pmat, int size)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			if( pmat[i*size+j].r != 0.0 || pmat[i*size+j].i != 0.0)
				return 1;
		}
	}

	return 0;
}

int ScanMatrixDiag(doublecomplex *pmat, int size)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			if(i!=j)
				if( pmat[i*size+j].r != 0.0 || pmat[i*size+j].i != 0.0)
					return 1;
			if(i==j)
				if(pmat[i*size+j].r != 2.0 || pmat[i*size+j].i != 0.0)
					return 1;
		}
	}

	return 0;
}
