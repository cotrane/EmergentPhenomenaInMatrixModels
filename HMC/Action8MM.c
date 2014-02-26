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
 * computation of action for d=8 model; both pure Yang-Mills-Myers model and modified YM-Myers model if \lambda is not 0
 */

double action8MM(doublecomplex *pmat[], double *trYM1, double *trCS1, double coupling, int size, int nummat, double C2, double C3, double lambda,
		double *Action0, double *Action1)
{
	int i,j, mat1=0, mat2, mat3, pos1, pos2;
	double mat_c[24];
	double trCS=0, trYM=0, trrand=0, trrand1=0, con;
	doublecomplex *pcomm[28];

	for(i=0;i<28;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	// computation of commutators and storing them in pcomm
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
	// squaring commutators
	for(i=0; i<28; i++)
	{
		diagMulti(&trYM, pcomm[i], pcomm[i], size, 1.0);
	}

	//Myers-Term
	// result saved in 'trCS'
	for(mat1=0; mat1<(nummat-1); mat1++)
	{
		// reading in values of structure constant that include matrix number 'mat1'; values saved in 'mat_c'
		thecode(mat1, mat_c);
		for(i=0;i<=18;i+=6)
		{
			if(mat_c[i]>=0.0)
			{
				mat2=mat_c[i];
				mat3=mat_c[i+1];
				con=mat_c[i+2];
				if(mat1<mat2 && mat2<mat3)
				{
					pos1=0;
					for(j=0;j<mat2;j++)
					{
						pos1+=((nummat-1)-j);
					}
					diagMulti(&trCS, pcomm[pos1+mat3-1-mat2], pmat[mat1], size, con);
				}
			}
			else
			{
				continue;
			}
		}
	}

	// multiplying terms by constants of model
	*trYM1 = size*(-1/2.0 *(trYM));
	*trCS1 = -size*coupling/pow(size,0.25)*(2.0 *(trCS));
	// adding terms up to final result in 'trrand' and for separate evaluation in 'Action0'
	trrand = *trYM1 + *trCS1;
	*Action0 = trrand;

	// if lambda is not 0 computation of additional term in modified YM-Myers model; saved in trrand1 and Action1 and added to trrand
	if(lambda!=0.0)
	{
		trrand1 += S1(pmat, size, C2, C3, lambda, coupling, nummat);
		*Action1 = trrand1;
		trrand += trrand1;
	}

	for(i=0;i<28;i++)
	{
		free(pcomm[i]);
	}

	return trrand;
}
