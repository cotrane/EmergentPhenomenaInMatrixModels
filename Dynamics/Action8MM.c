#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"
#include "MatrixMan.h"


double action8MM(doublecomplex *pmat[], double coupling, int size, int nummat, double C2, double C3, double lambda, double *Action0, double *Action1)
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

	for(i=0; i<28; i++)															//Squaring them
	{
		diagMulti(&trYM, pcomm[i], pcomm[i], size, 1.0);
	}

	//CS-Term

	for(mat1=0; mat1<(nummat-1); mat1++)
	{
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
	//Initial Energy

	trYM = size*(-1/2.0 *(trYM));
	trCS = -size*coupling/pow(size,0.25)*(2.0 *(trCS));

	trrand = trYM + trCS ;
	*Action0 = trrand;
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
