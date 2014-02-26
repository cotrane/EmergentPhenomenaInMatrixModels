/*
 * EVtoDist_Modulus_nH.c
 *
 *  Created on: 24 Jan 2013
 *      Author: tkaltenbrunner
 */

#include<stdio.h>
#include<stdlib.h>
#include <float.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "EV.h"

/*
 * computes EV distribution of modulus of non-Hermitian matrix; saves the EV distribution, not the individual EV's
 */
int EVtoDist_Modulus_nH(int **out, doublecomplex *in, int size, double binwidth, int *bins,
		double *LargeEVmod, double *SmallEVmod, int *numEVmod)
{
	int i, Bin, BinsEV, *tmp, numev2=0;
	static int init=0;
	static double *SavEVreal, *SavEVimag, *SavEVmod;
	double RangeEVtmp=0, LargeEVtmp2=0, SmallEVtmp2=0;
	doublecomplex LargeEVtmp, SmallEVtmp;

	if(init==0){
		SavEVreal = (double*) malloc(sizeof(double)*size);
		SavEVimag = (double*) malloc(sizeof(double)*size);
		SavEVmod = (double*) malloc(sizeof(double)*size);
		init=1;
	}
	memset(SavEVreal, 0, sizeof(double)*size);
	memset(SavEVimag, 0, sizeof(double)*size);
	memset(SavEVmod, 0, sizeof(double)*size);
	LargeEVtmp.r=0; LargeEVtmp.i=0; SmallEVtmp.r=0; SmallEVtmp.i=0;

	// compute EV's of non-Hermitian matrix using LAPACK routine
	SaveEV_nonHermitian(SavEVreal, SavEVimag, &LargeEVtmp, &SmallEVtmp, in, size, &numev2);

	/*
	 * computing distribution for rho(r) = w(r)/r with w(r) from \int dr r rho(r) = 1; r = \sqrt(x^2+y^2)
	 */

	// find largest and smallest EV
	for(i=0;i<size;i++)
	{
		SavEVmod[i] = sqrt( pow(SavEVreal[i],2) + pow(SavEVimag[i],2) );

		if( SavEVmod[i] < SmallEVtmp2)
		{
			SmallEVtmp2 = SavEVmod[i];
		}
		else if( SavEVmod[i] > LargeEVtmp2)
		{
			LargeEVtmp2 = SavEVmod[i];
		}
	}

	// find range of EV's
	if(*bins == 0)
	{
		RangeEVtmp = LargeEVtmp2 -SmallEVtmp2;
		BinsEV = RangeEVtmp/binwidth;

		*(out) = (int*) malloc(sizeof(int)*BinsEV);
		memset(*out, 0, sizeof(int)*BinsEV);

		*SmallEVmod = SmallEVtmp2;
		*LargeEVmod = LargeEVtmp2;
		*bins = BinsEV;
	}
	// extend upper range of previous EV dist if necessary
	if(LargeEVtmp2 > *LargeEVmod)
	{
		RangeEVtmp = LargeEVtmp2 - *SmallEVmod;
		BinsEV = RangeEVtmp/binwidth;

		tmp = (int*) malloc(sizeof(int)* *bins);
		memcpy(tmp, *out, sizeof(int)* *bins);
		free(*out);
		*out = (int*) malloc(sizeof(int)*BinsEV);
		memset(*out, 0, sizeof(int)*BinsEV);

		for(i=0;i< *bins;i++)
			(*out)[i] = tmp[i];
		free(tmp);

		*LargeEVmod = LargeEVtmp2;
		*bins = BinsEV;
	}
	// extend lower range of previous EV dist if necessary
	if(SmallEVtmp2 < *SmallEVmod)
	{
		RangeEVtmp = *LargeEVmod - SmallEVtmp2;
		BinsEV = RangeEVtmp/binwidth;

		tmp = (int*) malloc(sizeof(int)* *bins);
		memcpy(tmp, *out, sizeof(int)* *bins);
		free(*out);
		*out = (int*) malloc(sizeof(int)*BinsEV);
		memset(*out, 0, sizeof(int)*BinsEV);

		for(i=(BinsEV- *bins);i< BinsEV;i++)
			(*out)[i] = tmp[i-(BinsEV-*bins)];
		free(tmp);

		*SmallEVmod = SmallEVtmp2;
		*bins = BinsEV;
	}
	// update EV distribution
	for(i=0; i< size; i++)
	{
		Bin = (SavEVmod[i]- *SmallEVmod)/(binwidth*1.0);
		if(Bin > (*bins))
		{
			printf("Bin=%d\tbins=%d\tSavEVmod[%d]=%f\tSmallEVmod=%f\n", Bin, *bins, i, SavEVmod[i], *SmallEVmod);
			printf("too far!\n");
		}
		if(Bin==(*bins))
			Bin = *bins -1;
		(*out)[Bin] += 1;
	}

	*numEVmod = *numEVmod + size;

	return 0;
}


