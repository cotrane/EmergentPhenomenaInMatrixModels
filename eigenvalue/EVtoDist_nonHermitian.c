/*
 * EVtoDist_nonHermitian.c
 *
 *  Created on: 23 Jan 2013
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
 * computes EV distribution of real and imaginary part of non-Hermitian matrix; saves the EV distribution, not the individual EV's
 */
int EVtoDist_nonHermitian(int **out1, int **out2, doublecomplex *in, int size, double binwidth, int bins[2],
		doublecomplex *LargeEV,	doublecomplex *SmallEV, int *numev)
{
	int i, Bin, BinsEV[2], *tmp, numev2=0;
	static int init=0;
	static double *SavEVreal, *SavEVimag;
	doublecomplex LargeEVtmp, SmallEVtmp, RangeEVtmp;

	if(init==0){
		SavEVimag = (double*) malloc(sizeof(double)*size);
		SavEVreal = (double*) malloc(sizeof(double)*size);
		init=1;
	}
	memset(SavEVreal, 0, sizeof(double)*size);
	memset(SavEVimag, 0, sizeof(double)*size);
	LargeEVtmp.r=0; LargeEVtmp.i=0; SmallEVtmp.r=0; SmallEVtmp.i=0;

	// compute EV's of non-Hermitian matrix using LAPACK routine
	SaveEV_nonHermitian(SavEVreal, SavEVimag, &LargeEVtmp, &SmallEVtmp, in, size, &numev2);

	/*
	 * computing distribution of real part
	 */

	// find largest and smallest EV
	// find range of EV's
	if(bins[0] == 0)
	{
		RangeEVtmp.r = LargeEVtmp.r-SmallEVtmp.r;
		BinsEV[0] = RangeEVtmp.r/binwidth;

		*(out1) = (int*) malloc(sizeof(int)*BinsEV[0]);
		memset(*out1, 0, sizeof(int)*BinsEV[0]);

		SmallEV->r = SmallEVtmp.r;
		LargeEV->r = LargeEVtmp.r;
		bins[0] = BinsEV[0];
	}
	// extend upper range of previous EV dist if necessary
	if(LargeEVtmp.r > (LargeEV->r))
	{
		RangeEVtmp.r = LargeEVtmp.r - SmallEV->r;
		BinsEV[0] = RangeEVtmp.r/binwidth;

		tmp = (int*) malloc(sizeof(int)* bins[0]);
		memcpy(tmp, *out1, sizeof(int)* bins[0]);
		free(*out1);
		*out1 = (int*) malloc(sizeof(int)*BinsEV[0]);
		memset(*out1, 0, sizeof(int)*BinsEV[0]);

		for(i=0;i< bins[0];i++)
			(*out1)[i] = tmp[i];
		free(tmp);

		LargeEV->r = LargeEVtmp.r;
		bins[0] = BinsEV[0];
	}
	// extend lower range of previous EV dist if necessary
	if(SmallEVtmp.r < (SmallEV->r))
	{
		RangeEVtmp.r = LargeEV->r - SmallEVtmp.r;
		BinsEV[0] = RangeEVtmp.r/binwidth;

		tmp = (int*) malloc(sizeof(int)* bins[0]);
		memcpy(tmp, *out1, sizeof(int)* bins[0]);
		free(*out1);
		*out1 = (int*) malloc(sizeof(int)*BinsEV[0]);
		memset(*out1, 0, sizeof(int)*BinsEV[0]);

		for(i=0;i< bins[0];i++)
			(*out1)[BinsEV[0]-i-1] = tmp[bins[0]-i-1];
		free(tmp);

		SmallEV->r = SmallEVtmp.r;
		bins[0] = BinsEV[0];
	}
	// update EV distribution
	for(i=0; i< size; i++)
	{
		Bin = (SavEVreal[i]- SmallEV->r)/(binwidth*1.0);
		if(Bin>bins[0])
		{
			printf("Bin=%d\tbins[0]=%d\tSavEVreal[%d]=%f\tSmallEV.r=%f\tLargeEV.r=%f\n", Bin, bins[0], i, SavEVreal[i], SmallEV->r, LargeEV->r);
			printf("too far!\n");
		}
		if(Bin==(bins[0]))
			Bin = bins[0]-1;
		(*out1)[Bin] += 1;
	}

	/*
	 * computing distribution of imaginary part
	 */

	// find largest and smallest EV
	// find range of EV's
	if(bins[1] == 0)
	{
		RangeEVtmp.i = LargeEVtmp.i-SmallEVtmp.i;
		BinsEV[1] = RangeEVtmp.i/binwidth;

		*out2 = (int*) malloc(sizeof(int)*BinsEV[1]);
		memset(*out2, 0, sizeof(int)*BinsEV[1]);

		SmallEV->i = SmallEVtmp.i;
		LargeEV->i = LargeEVtmp.i;
		bins[1] = BinsEV[1];
	}
	// extend upper range of previous EV dist if necessary
	if(LargeEVtmp.i > LargeEV->i)
	{
		RangeEVtmp.i = LargeEVtmp.i - SmallEV->i;
		BinsEV[1] = RangeEVtmp.i/binwidth;

		tmp = (int*) malloc(sizeof(int)* bins[1]);
		memcpy(tmp, *out2, sizeof(int)* bins[1]);
		free(*out2);
		*out2 = (int*) malloc(sizeof(int)*BinsEV[1]);
		memset(*out2, 0, sizeof(int)*BinsEV[1]);

		for(i=0;i< bins[1];i++)
			(*out2)[i] = tmp[i];
		free(tmp);

		LargeEV->i = LargeEVtmp.i;
		bins[1] = BinsEV[1];
	}
	// extend lower range of previous EV dist if necessary
	if(SmallEVtmp.i < SmallEV->i)
	{
		RangeEVtmp.i = LargeEV->i - SmallEVtmp.i;
		BinsEV[1] = RangeEVtmp.i/binwidth;

		tmp = (int*) malloc(sizeof(int)* bins[1]);
		memcpy(tmp, *out2, sizeof(int)* bins[1]);
		free(*out2);
		*out2 = (int*) malloc(sizeof(int)*BinsEV[1]);
		memset(*out2, 0, sizeof(int)*BinsEV[1]);

		for(i=0;i< bins[1];i++)
			(*out2)[BinsEV[1]-i-1] = tmp[bins[1]-i-1];
		free(tmp);

		SmallEV->i = SmallEVtmp.i;
		bins[1] = BinsEV[1];
	}
	// update EV distribution
	for(i=0; i< size; i++)
	{
		Bin = (SavEVimag[i]- SmallEV->i)/(binwidth*1.0);
		if(Bin>bins[1])
		{
			printf("Bin=%d\tbins[1]=%d\tSavEVimag[%d]=%f\tSmallEV.i=%f\n", Bin, bins[1], i, SavEVimag[i], SmallEV->i);
			printf("too far!\n");
		}
		if(Bin==(bins[1]))
			Bin = bins[1]-1;
		(*out2)[Bin] += 1;
	}

	*numev = *numev+size;

	return 0;
}


