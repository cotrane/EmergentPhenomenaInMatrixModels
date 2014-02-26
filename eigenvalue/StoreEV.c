#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

extern char **pfilenamest;
extern FILE **poutst;
char pfilenamelambdast[256];
FILE *poutlambdast;

void StoreEV(double *pSavEV[], int NUMMAT, int numev[])
{
	int i,j;

	if(NUMMAT!=1)
	{
		for(i=0;i<NUMMAT;i++)
		{
			for(j=0;j<numev[i];j++)
			{
				fprintf(poutst[i], "%lf\n", (*(pSavEV[i]+j)));
			}
			fflush(poutst[i]);
		}
		for(i=0;i<NUMMAT;i++)
		{
			memset(pSavEV[i],0,sizeof(double)*numev[i]);
			numev[i]=0;
		}
	}

	if(NUMMAT==1)
	{
		i=0;
		for(j=0;j<numev[i];j++)
		{
			fprintf(poutlambdast, "%lf\n", (*(pSavEV[i]+j)));
		}
		fflush(poutlambdast);
		for(i=0;i<NUMMAT;i++)
		{
			memset(pSavEV[i],0,sizeof(double)*numev[i]);
			numev[i]=0;
		}
	}
}

