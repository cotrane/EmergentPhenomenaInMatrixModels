#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"


int dtensor(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE)
{
	int i,j,k;
	double tr=0;
	doublecomplex *AntiCommm[1];
	FILE *out;

	out = fopen("/home/tkaltenbrunner/workspace/symmtensor.txt", "w");
	AntiCommm[0] = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);

	for(i=0;i<NUMMAT;i++)
	{
		for(j=0;j<NUMMAT;j++)
		{
			for(k=0;k<NUMMAT;k++)
			{
				memset(AntiCommm[0], 0, MATRIX_SIZE*MATRIX_SIZE*sizeof(doublecomplex));
				tr=0;
				AntiComm(AntiCommm[0], pmat[i], pmat[j], MATRIX_SIZE, 1.0);
				diagMulti(&tr, AntiCommm[0], pmat[k], MATRIX_SIZE, 2.0);
				if(tr<-0.001 || tr>0.001)
				{
					printf("%d,%d,%d : %f\n", i,j,k, tr);
					fprintf(out, "%f\n", tr);
				}
			}
		}
	}

	fclose(out);
	return 0;
}
