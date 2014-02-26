#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "EV.h"
#include <float.h>

/*
 * computes the eigenvalues of a matrix and saves its distribution - not the individual eigenvalues
 */

int EVonthefly(int **out, doublecomplex *in, int size, double binwidth, int *bins, double *LargeEV, double *SmallEV, int *numev)
{
	int i, Bin, BinsEV, *tmp, numev2=0;
	double *SavEV, LargeEVtmp=0, SmallEVtmp=0, RangeEVtmp=0;

	SavEV = (double*) malloc(sizeof(double)*size);
	memset(SavEV, 0, sizeof(double)*size);

	// compute eigenvalues of matrix using EV routine from LAPACK
	SaveEV(SavEV, &LargeEVtmp, &SmallEVtmp, in, size, &numev2);

	// compute range of EV's, largest and smalles EV
	if(*bins == 0)
	{
		RangeEVtmp = LargeEVtmp-SmallEVtmp;
		BinsEV = RangeEVtmp/binwidth;

		*out = (int*) malloc(sizeof(int)*BinsEV);
		memset(*out, 0, sizeof(int)*BinsEV);

		*SmallEV = SmallEVtmp;
		*LargeEV = LargeEVtmp;
		*bins = BinsEV;
	}

	// extend range of array for EV distribution if largest current EV is larger than previous largest EV
	if(LargeEVtmp > *LargeEV)
	{
		RangeEVtmp = LargeEVtmp - *SmallEV;
		BinsEV = RangeEVtmp/binwidth;

		tmp = (int*) malloc(sizeof(int)* *bins);
		memcpy(tmp, *out, sizeof(int)* *bins);
		free(*out);
		*out = (int*) malloc(sizeof(int)*BinsEV);
		memset(*out, 0, sizeof(int)*BinsEV);

		for(i=0;i< *bins;i++)
			(*out)[i] = tmp[i];
		free(tmp);

		*LargeEV = LargeEVtmp;
		*bins = BinsEV;
	}
	// extend range of array for EV distribution if smallest current EV is smaller than previous smallest EV
	if(SmallEVtmp < *SmallEV)
	{
		RangeEVtmp = *LargeEV - SmallEVtmp;
		BinsEV = RangeEVtmp/binwidth;

		tmp = (int*) malloc(sizeof(int)* *bins);
		memcpy(tmp, *out, sizeof(int)* *bins);
		free(*out);
		*out = (int*) malloc(sizeof(int)*BinsEV);
		memset(*out, 0, sizeof(int)*BinsEV);

		for(i=(BinsEV- *bins);i< BinsEV;i++)
			(*out)[i] = tmp[i-(BinsEV-*bins)];
		free(tmp);

		*SmallEV = SmallEVtmp;
		*bins = BinsEV;
	}

	// generate distribution
	for(i=0; i< size; i++)
	{
		Bin = (SavEV[i]- *SmallEV)/(binwidth*1.0);
		if(Bin>*bins)
		{
			printf("Bin=%d\tbins=%d\tSavEV[%d]=%f\tSmallEV=%f\n", Bin, *bins, i, SavEV[i], *SmallEV);
			printf("too far!\n");
		}
		if(Bin==(*bins))
			Bin = *bins-1;
		(*out)[Bin] += 1;
	}

	*numev = *numev+size;
	free(SavEV);

	return 0;
}



