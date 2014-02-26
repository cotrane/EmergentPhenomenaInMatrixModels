#include<f2c.h>
#include "prng.h"
#include"hmc.h"
#include "RandomGens.h"

int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE)
{
	int i,m, n;
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
					(pmat[i]+(m*MATRIX_SIZE + n))->r = prng_get_double();
					lastel += (pmat[i]+(m*MATRIX_SIZE + n))->r;
				}
				else
				{
					(pmat[i]+(m*MATRIX_SIZE + n))->r = prng_get_double();
					(pmat[i]+(m*MATRIX_SIZE + n))->i = prng_get_double();

					(pmat[i]+(n*MATRIX_SIZE + m))->r = (pmat[i]+(m*MATRIX_SIZE + n))->r;
					(pmat[i]+(n*MATRIX_SIZE + m))->i = -(pmat[i]+(m*MATRIX_SIZE + n))->i;
				}
			}
		}
		(pmat[i]+(MATRIX_SIZE-1)*MATRIX_SIZE + (MATRIX_SIZE-1))->r = -lastel;
		(pmat[i]+(MATRIX_SIZE-1)*MATRIX_SIZE + (MATRIX_SIZE-1))->i = 0.0;

	}

	return 0;
}

int gen_randphicplxtrace(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE)
{
	int i,m, n;

	for(i=0; i<NUMMAT; i++)
	{
		for(m = 0; m < (MATRIX_SIZE) ;m++)
		{
			for(n=m; n < MATRIX_SIZE ; n++)
			{
				if(m == n)
				{
					(pmat[i]+(m*MATRIX_SIZE + n))->r = prng_get_double();
				}
				else
				{
					(pmat[i]+(m*MATRIX_SIZE + n))->r = prng_get_double();
					(pmat[i]+(m*MATRIX_SIZE + n))->i = prng_get_double();

					(pmat[i]+(n*MATRIX_SIZE + m))->r = (pmat[i]+(m*MATRIX_SIZE + n))->r;
					(pmat[i]+(n*MATRIX_SIZE + m))->i = -(pmat[i]+(m*MATRIX_SIZE + n))->i;
				}
			}
		}

	}

	return 0;
}

int gen_randphicplxalt(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE)
{
	int i,j,m, n;
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
					(pmat[i]+(m*MATRIX_SIZE + n))->r = prng_get_double();
					lastel += (pmat[i]+(m*MATRIX_SIZE + n))->r;
				}
				else
				{
					(pmat[i]+(m*MATRIX_SIZE + n))->r = prng_get_double();
					(pmat[i]+(m*MATRIX_SIZE + n))->i = prng_get_double();

					(pmat[i]+(n*MATRIX_SIZE + m))->r = (pmat[i]+(m*MATRIX_SIZE + n))->r;
					(pmat[i]+(n*MATRIX_SIZE + m))->i = -(pmat[i]+(m*MATRIX_SIZE + n))->i;
				}
			}
		}
		for(j=0;j<MATRIX_SIZE;j++)
		{
			(pmat[i]+(j*MATRIX_SIZE + j))->r -= lastel/MATRIX_SIZE;
		}
	}

	return 0;
}

int gen_randphicplxdiagtrace(doublecomplex pmat[], int MATRIX_SIZE, double mult)
{
	int m;

	for(m = 0; m < MATRIX_SIZE ; m++)
	{
		pmat[m*MATRIX_SIZE+m].r = mult*gauss_randomnr(1);
	}

	return 0;
}
