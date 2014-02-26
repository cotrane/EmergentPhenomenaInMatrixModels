#include<f2c.h>
#include "RandomGens.h"

int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE)
{
	int i,m, n;
	double lastel;

	for(i=0; i<NUMMAT; i++)
	{
		lastel=0.0;
		for(m = 0; m < MATRIX_SIZE ;m++)
		{
			for(n=m; n < MATRIX_SIZE ; n++)
			{
				if(m == n && n<(MATRIX_SIZE-1))
				{
					(pmat[i]+(m*MATRIX_SIZE + n))->r = mt_ldrand();
					lastel += (pmat[i]+(m*MATRIX_SIZE + n))->r;
				}
				else
				{
					(pmat[i]+(m*MATRIX_SIZE + n))->r = mt_ldrand();
					(pmat[i]+(m*MATRIX_SIZE + n))->i = mt_ldrand();

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
