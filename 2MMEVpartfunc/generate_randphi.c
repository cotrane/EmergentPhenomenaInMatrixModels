#include<f2c.h>
#include<math.h>
#include "RandomGens.h"

int gen_randphicplx(double pmat[], int MATRIX_SIZE)
{
	int m;
	double lastel=0.0;

	for(m = 0; m < (MATRIX_SIZE-1) ; m++)
	{
		pmat[m] = mt_ldrand();
		lastel += pmat[m];
	}
	pmat[MATRIX_SIZE-1] = -lastel;

	return 0;
}

int gen_randphicplxtrace(double pmat[], int MATRIX_SIZE, double mult)
{
	int m;

	for(m = 0; m < MATRIX_SIZE ; m++)
	{
		pmat[m] = mult*mt_ldrand()+1.0;
	}

	return 0;
}

int gen_gaussMat_cplx_trace(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a)
{
	int i, m, n;

	for(i=0; i<NUMMAT; i++)
	{
		for(m = 0; m < (MATRIX_SIZE) ;m++)
		{
			for(n=m; n < MATRIX_SIZE ; n++)
			{
				if(m == n)
				{
					(pmom[i]+(m*MATRIX_SIZE + n))->r = a*gauss_randomnr(1);
					(pmom[i]+(m*MATRIX_SIZE + n))->i = 0.0;
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
	}
	return 0;
}
