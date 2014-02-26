#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

/*
 * compute EV distribution from EV's saved in pSavEV; those were computed with SaveEV
 */
extern FILE *pout;

void DataDistribution(double *pSavEV, double LargEV, double SmallEV, double BINWIDTH, int num, int numev)
{
	int i;
	double RangeEV, interval;
	int *pev, Bin, BinsEV;

	RangeEV = LargEV - SmallEV;
	BinsEV = RangeEV/BINWIDTH +1;
	pev = (int*) malloc(BinsEV * sizeof(int));
	memset(pev, 0, BinsEV * sizeof(int));

	for(i=0; i< numev; i++)
	{
		Bin = (pSavEV[i]- SmallEV)/(BINWIDTH*1.0);
		if(Bin>BinsEV)
		{
			printf("too far!\n");
		}
		if(Bin==(BinsEV))
			Bin = BinsEV-1;
		pev[Bin] += 1;
	}

	for(i=0; i<BinsEV;i++)
	{
		interval = SmallEV + BINWIDTH*i;
		fprintf(pout, "%lf\t%lf\n", interval, pev[i]/(numev*BINWIDTH*1.0));
		fflush(pout);
	}

	free(pev); // he wants to break free
}

