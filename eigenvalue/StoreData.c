#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

/*
 * stores eigenvalues saved in pSavEV into file opened with FileOpen.c
 */
extern FILE *poutst;

void StoreData(double *pSavEV, int *numev)
{
	int j;

	for(j=0;j<*numev;j++)
	{
		fprintf(poutst, "%lf\n", pSavEV[j]);
	}
	fflush(poutst);

	memset(pSavEV,0,sizeof(double)* *numev);
	*numev=0;
}

