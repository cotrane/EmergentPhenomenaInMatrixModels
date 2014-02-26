#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "EV.h"

/*
 * compute EV's of matrix and store those EV's in array pSavEV
 */
void SaveEV(double *pSavEV, double *LargEV, double *SmallEV, doublecomplex *pmat, int MATRIX_SIZE, int *numev)
{
	int i;
	int N=MATRIX_SIZE*MATRIX_SIZE;
	doublecomplex *pmatev;
	doublereal *pvec;

	pvec = (doublereal*) malloc(sizeof(doublecomplex)*MATRIX_SIZE);
	pmatev = (doublecomplex*) malloc(sizeof(doublecomplex)*N);

	memcpy(pmatev, pmat, N*sizeof(doublecomplex));
	findEigenvalues(pvec, pmatev, MATRIX_SIZE);

	for(i=0; i<MATRIX_SIZE; i++)
	{
		pSavEV[*numev] = pvec[i];
		(*numev)++;

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
/*
 * compute EV's of non-Hermitian matrix and store those EV's in array pSavEVr and pSavEVi (real and imaginary part respectively)
 */
void SaveEV_nonHermitian(double *pSavEVr, double *pSavEVi, doublecomplex *LargEV, doublecomplex *SmallEV, doublecomplex *pmat,
		int MATRIX_SIZE, int *numev)
{
	int i;
	int N=MATRIX_SIZE*MATRIX_SIZE;
	static int init=0;
	static doublecomplex *pmatev, *pvec;

	if(init==0){
		pvec = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE);
		pmatev = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		init=1;
	}
	memcpy(pmatev, pmat, N*sizeof(doublecomplex));
	findEigenvalues_NonHermitian(pvec, pmatev, MATRIX_SIZE);

	for(i=0; i<MATRIX_SIZE; i++)
	{
		pSavEVr[*numev] = pvec[i].r;
		pSavEVi[*numev] = pvec[i].i;
		(*numev)++;

		if( (pvec[i].r) < SmallEV->r)
		{
			SmallEV->r = pvec[i].r;
		}
		else if( (pvec[i].r) > LargEV->r)
		{
			LargEV->r = pvec[i].r;
		}

		if( (pvec[i].i) < SmallEV->i)
		{
			SmallEV->i = pvec[i].i;
		}
		else if( (pvec[i].i) > LargEV->i)
		{
			LargEV->i = pvec[i].i;
		}
	}

}
