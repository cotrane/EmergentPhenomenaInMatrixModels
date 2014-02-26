#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include <time.h>

double computeAutocorrelationtime(double *Data, double *Mean, int loops, double *Sigmast2);
int computeJackknife(double *SigmaJK, double *SigmaJKspecheat, double *Data, double loops, double *Mean, double totEnergy, double *Tau);

float Errors(double *pMean, double *pSigmast, double *pSigmast2, double *pSigmaJK, double *pSigmaJKspecheat, double *pTau, double *pEnergy, int size, int kint)
{
	int i;
	double k = (kint)*1.0;
	double Mean2 = 0;
	double EnergyJK[kint];
	double totEnergy=0;

	memset(EnergyJK, 0, kint*sizeof(double));

	for(i=0; i<kint; i++)				//calculates mean value, square of it and the sum of the energies of all timesteps for use in jackknife
	{
		*pMean += pEnergy[i]/k;
		Mean2 += (pEnergy[i]*pEnergy[i])/k;
		totEnergy += pEnergy[i];
	}

	for(i=0;i<kint;i++)				//computes square of standard deviation and square of difference between energy of diff steps and mean value for autocorrelation
	{
		*pSigmast2 += (pEnergy[i] - *pMean)*(pEnergy[i] - *pMean)/k;
	}
	*pSigmast = sqrt(*pSigmast2/(k-1.0));		//standard deviation without any correction coefficient

	//Autocorrelation

	*pTau = computeAutocorrelationtime(pEnergy, pMean, kint, pSigmast2);

	//Jackknife

	computeJackknife(pSigmaJK, pSigmaJKspecheat, pEnergy, kint, pMean, totEnergy, pTau);

	printf("SpecHeat = %lf\tSigmaJKsh = %lf\n", *pSigmast2/(size*size), *pSigmaJKspecheat);
	printf("Mean = %lf\tSigmast = %lf\tSigmaJK = %lf\tTau = %lf\n",	*pMean, *pSigmast, *pSigmaJK, *pTau);

	return 0;
}

int computeJackknife(double *SigmaJK, double *SigmaJKspecheat, double *Data, double loops, double *Mean, double totEnergy, double *Tau)
{
	int d, i, j, numd;
	double SigmaJKtry, *EnergyJK, *EnergyJK2, SigmaJKspecheattry, *Specheattry, meanSpecheat=0;

	EnergyJK = (double*) malloc(sizeof(double)*(loops+1));
	EnergyJK2 = (double*) malloc(sizeof(double)*(loops+1));
	Specheattry = (double*) malloc(sizeof(double)*(loops+1));
	//Jackknife

	for(d=0; d<500; d+=5)			//how many elements are omitted
	{
		while(d<0)
		{
			d=d+5;
		}
		memset(Specheattry, 0, (loops+1)*sizeof(double));
		SigmaJKtry = 0.0;
		SigmaJKspecheattry = 0.0;
		meanSpecheat = 0.0;

		numd = loops/d;
		for(j=0; j<numd; j++)
		{
			memset(EnergyJK, 0, (loops+1)*sizeof(double));
			memset(EnergyJK2, 0, (loops+1)*sizeof(double));
			for(i=0; i<d; i++)			//summing over omitted elements
			{
				EnergyJK[j] += Data[d*j+i];
				EnergyJK2[j] += Data[d*j+i]*Data[d*j+i];
			}

			EnergyJK[j] = (totEnergy - EnergyJK[j]) /(loops-d*1.0);									//subtracting omitted elements from total sum of energies
			EnergyJK2[j] = (totEnergy*totEnergy - EnergyJK2[j]) /(loops-d*1.0);						//subtracting square of omitted elements from tot sum of energies for spec.heat
			SigmaJKtry +=  (numd - 1.0)*(EnergyJK[j] - *Mean)*(EnergyJK[j] - *Mean) /numd;			//computing corrected square of sigma

			Specheattry[j] = (EnergyJK2[j] - EnergyJK[j]*EnergyJK[j])/numd;							//computing specheat without omitted elements
			meanSpecheat += Specheattry[j]/numd;													//computing mean specheat
		}
		for(j=0;j<numd;j++)
		{
			SigmaJKspecheattry += (numd - 1.0)*(Specheattry[j]-meanSpecheat)*(Specheattry[j]-meanSpecheat)/numd;
		}
		SigmaJKtry = sqrt(SigmaJKtry);			//Jackknife corrected standard deviation
		if(SigmaJKtry > *SigmaJK)
		{
			*SigmaJK = SigmaJKtry;				//saving largest value of sigma corresponding to specific # of omitted elements
			*SigmaJKspecheat = sqrt(SigmaJKspecheattry);				//Jackknife corrected standard deviation for specheat
		}
	}

	return 0;
}

double computeAutocorrelationtime(double *Data, double *Mean, int loops, double *Sigmast2)
{
	int i,d, border;
	double Autoc, Tau = 0.5;

	border=loops/16;
	if(border>5000)
	{
		border=5000;
	}
	for(d=1; d<border; d++)
	{
		Autoc=0;
		for(i=0; i<(loops-d); i++)
		{
			Autoc += (Data[i]-*Mean)*(Data[i+d]-*Mean)/(loops*1.0/*-(d*1.0)*/);
		}
		Autoc /= *Sigmast2;
		Tau += Autoc;
	}

	return Tau;
}
