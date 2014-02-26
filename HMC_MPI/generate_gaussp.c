#include<f2c.h>
#include<stdio.h>
#include<math.h>
#include "hmc.h"
#include "montecarlo.h"
#include "RandomGens.h"


int gen_gaussmomcplx(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a)
{
	int i, m, n;
	double lastel;

	for(i=0; i<NUMMAT; i++)
	{
		lastel=0.0;
		for(m = 0; m < (MATRIX_SIZE-1) ;m++)
		{
			for(n=m; n < MATRIX_SIZE ; n++)
			{
				if(m == n)
				{
					(pmom[i]+(m*MATRIX_SIZE + n))->r = a*gauss_randomnr(1);
					(pmom[i]+(m*MATRIX_SIZE + n))->i = 0.0;
					lastel += (pmom[i]+(m*MATRIX_SIZE + n))->r;
				}
				else
				{
 					(pmom[i]+(m*MATRIX_SIZE + n))->r = a/(sqrt(2.0)) * gauss_randomnr(1);
					(pmom[i]+(m*MATRIX_SIZE + n))->i = a/(sqrt(2.0)) * gauss_randomnr(1);
					(pmom[i]+(n*MATRIX_SIZE + m))->r = (pmom[i]+(m*MATRIX_SIZE + n))->r;
					(pmom[i]+(n*MATRIX_SIZE + m))->i = -(pmom[i]+(m*MATRIX_SIZE + n))->i;
				}
			}
		}
		(pmom[i]+(MATRIX_SIZE-1)*MATRIX_SIZE + (MATRIX_SIZE-1))->r = -lastel;
		(pmom[i]+(MATRIX_SIZE-1)*MATRIX_SIZE + (MATRIX_SIZE-1))->i = 0.0;
	}
	return 0;
}

int gen_gaussmomcplxalt(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a)
{
	int i, j, m, n;
	double lastel;

	for(i=0; i<NUMMAT; i++)
	{
		lastel=0.0;
		for(m = 0; m < MATRIX_SIZE ;m++)
		{
			for(n=m; n < MATRIX_SIZE ; n++)
			{
				if(m == n)
				{
					(pmom[i]+(m*MATRIX_SIZE + n))->r = a*gauss_randomnr(1);
					(pmom[i]+(m*MATRIX_SIZE + n))->i = 0.0;
					lastel += (pmom[i]+(m*MATRIX_SIZE + n))->r;
				}
				else
				{
 					(pmom[i]+(m*MATRIX_SIZE + n))->r = a/(sqrt(2.0)) * gauss_randomnr(1);
					(pmom[i]+(m*MATRIX_SIZE + n))->i = a/(sqrt(2.0)) * gauss_randomnr(1);
					(pmom[i]+(n*MATRIX_SIZE + m))->r = (pmom[i]+(m*MATRIX_SIZE + n))->r;
					(pmom[i]+(n*MATRIX_SIZE + m))->i = -(pmom[i]+(m*MATRIX_SIZE + n))->i;
				}
			}
		}
		for(j=0;j<MATRIX_SIZE;j++)
		{
			(pmom[i]+(j*MATRIX_SIZE + j))->r -= lastel/(MATRIX_SIZE*1.0);
		}
	}
	return 0;
}
