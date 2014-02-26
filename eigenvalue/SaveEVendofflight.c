#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

extern FILE *pout;

/*
 * together with EVonthefly; saves the EV distribution in a file
 */
int SaveEVendofflight(int *in, int BinsEV, double binwidth, double SmallEV, int numev)
{
	int i;
	double interval;

	for(i=0; i<BinsEV;i++)
	{
		interval = SmallEV + binwidth*i;
		fprintf(pout, "%lf\t%lf\n", interval, in[i]/(numev*binwidth*1.0));
		fflush(pout);
	}

	return 0;
}

