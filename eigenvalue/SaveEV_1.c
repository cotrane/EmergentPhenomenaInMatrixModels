#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "EV.h"

/*
 * computes eigenvalues and eigenvectors of Hermitian matrix;
 * EV's are saved in pSavEV,
 * eigenvectors saved in 'eigenvectors'
 * largest and smallest EV's saved in LargEV and SmallEV
 */
void SaveEV_1vector(double pSavEV[], double *LargEV, double *SmallEV, doublecomplex *eigenvectors, int MATRIX_SIZE)
{
	int i;
	int N=MATRIX_SIZE*MATRIX_SIZE;
	doublecomplex *pmatev;
	doublereal *pvec;

	pvec = (doublereal*) malloc(sizeof(doublecomplex)*MATRIX_SIZE);
	pmatev = (doublecomplex*) malloc(sizeof(doublecomplex)*N);

	memcpy(pmatev, eigenvectors, sizeof(doublecomplex)*N);
	findEigenvalues_vectors(pvec, pmatev, MATRIX_SIZE);
	memcpy(eigenvectors, pmatev, N*sizeof(doublecomplex));

	for(i=0; i<MATRIX_SIZE; i++)
	{
		pSavEV[i] = pvec[i];

		if((pvec[i])<*SmallEV)
		{
			*SmallEV = pvec[i];
		}
		else if(pvec[i] > *LargEV)
		{
			*LargEV = pvec[i];
		}
	}

	free(pvec);
	free(pmatev);
}
