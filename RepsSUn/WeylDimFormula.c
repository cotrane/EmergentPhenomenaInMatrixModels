/*
 * WeylDimFormula.c
 *
 *  Created on: 14 Nov 2012
 *      Author: tkaltenbrunner
 */

#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>

/*
 * computes the dimensionality of the matrix in terms of output is 'L',
 * representing boxes in Young-diagram; derived from Matrix Size 'size'
 * Weyls Dimension Formula valid for lie groups A_n, D_n, E_n
 */
int WeylDimFormula(int size, int liegroup)
{
	int j,k,l;
	int dim=0, dim1=0, sum, *rep;

	rep = (int*) malloc(sizeof(int)*(liegroup-1));
	memset(rep, 0, sizeof(int)*(liegroup-1));

	while(dim<=size)
	{
		dim1=dim;
		dim=1;
		sum=0;
		rep[0]+=1;
		for(j=1;j<liegroup;j++)
		{
			for(k=0;k<(liegroup-j);k++)
			{
				sum=0;
				for(l=k;l<(j+k);l++)
				{
					sum += rep[l];
				}
				dim = dim*(sum+j)/(1.*j);
			}
		}
	}
	rep[0]-=1;

	return rep[0];
}
