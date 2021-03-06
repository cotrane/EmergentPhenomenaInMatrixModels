#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

/*
 * routines to compute Error estimations:
 *    -) Jackknife Method
 *    -) integrated Autocorrelation time
 * parameters are either read from command line via ReadInputCommandLine or from file ErrorEval.txt via ReadInputFile.c
 */
int size, position, loops;

int computeJackknife(double *SigmaJK, double *SigmaJKspecheat, double *Data, double loops, double Mean, double totEnergy);
double computeAutocorrelationtime(double *Data, double Mean, int loops, double Sigmast2);
int ReadInputCL(int argc, char *argv[], int *step, int *size, double *alpha, int *loops, int *pos, char *pathout, char *pathin);
int ReadInputFile(int *step, int *size, double *alpha, int *loops, int *position, char *pathout, char *pathin);

int main(int argc, char *argv[])
{
	int i,count=0, countlines=1;
	char stream[30], pathin[1000], pathout[1000];
	FILE *in, *out;
	char filenamein[256], filenameout[256];
	double alpha, *Energy, Mean=0, Sigmast=0, Sigmast2=0, SigmaJK=0, SigmaJKspecheat, Tau=0;
	int step=1, lines=0;
	double k;
	double totEnergy=0;
	double specheat=0;

	// read in input parameter either from command line or from file
	if(argc>1)
		ReadInputCL(argc, argv, &step, &size, &alpha, &loops, &position, pathout, pathin);
	else
		ReadInputFile(&step, &size, &alpha, &loops, &position, pathout, pathin);

	printf("step=%d\tsize=%d\tloops=%d\tposition=%d\tpathout=%s\tpathin=%s\n", step, size, loops, position, pathout, pathin);

	sprintf(filenamein, "%s", pathin);
	sprintf(filenameout, "%s", pathout);
	in = fopen(filenamein,"r");
	out = fopen(filenameout, "a");

	Energy = (double*) malloc(sizeof(double)*(loops+1));
	memset(Energy, 0, sizeof(double)*(loops+1));

//	save the part of the input data that shall be used for error computation in array 'Energy'
	fseek(in, 0, SEEK_SET);
	while(fgets(stream, 30, in)!=NULL)
	{
		if(countlines == position)
		{
			while(lines!=loops)
			{
				if(lines%step==0)
				{
					Energy[count] = atof(stream);
					count++;
				}
				lines++;
				fgets(stream, 30, in);
			}
			break;
		}
		countlines++;
	}
	//setting loops to loops/steps
	loops = count;
	k = (loops)*1.0;

	//calculates mean value, square of it and the sum of the energies of all timesteps for use in jackknife
	for(i=0; i<loops; i++)
	{
		Mean += Energy[i]/k;
		totEnergy += Energy[i];
	}

	//computes square of standard deviation and square of difference between energy of diff steps and mean value for autocorrelation
	for(i=0;i<loops;i++)
	{
		Sigmast2 += (Energy[i] - Mean)*(Energy[i] - Mean)/(k);
	}
	//standard deviation without any correction coefficient
	Sigmast = sqrt(Sigmast2/(k-1));

	//Jackknife
	computeJackknife(&SigmaJK, &SigmaJKspecheat, Energy, loops, Mean, totEnergy);

	//Autocorrelation
	Tau = computeAutocorrelationtime(Energy, Mean, loops, Sigmast2);

	specheat = Sigmast2/(size*size);

	printf("spec. heat = %lf\tsigma spec. heat = %lf\n", specheat, SigmaJKspecheat);
	printf("Mean = %lf\tSigmast = %lf\tSigmaJK = %lf\tTau = %lf\n",	Mean, Sigmast, SigmaJK, Tau);
	fprintf(out, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", alpha, Mean, SigmaJK, specheat, SigmaJKspecheat, Tau);

	fclose(out);
	fclose(in);

	return 0;

}

/*
 * computes error estimate corrected for Autocorrelation using Jackknife Method
 */
int computeJackknife(double *SigmaJK, double *SigmaJKspecheat, double *Data, double loops, double Mean, double totEnergy)
{
	int d, i, j, numd;
	double SigmaJKtry, *EnergyJK, *EnergyJK2, SigmaJKspecheattry, *Specheattry, meanSpecheat=0;

	EnergyJK = (double*) malloc(sizeof(double)*(loops+1));
	EnergyJK2 = (double*) malloc(sizeof(double)*(loops+1));
	Specheattry = (double*) malloc(sizeof(double)*(loops+1));

	//Jackknife

	//how many elements are omitted
	for(d=1; d<500; d++)
	{
		memset(Specheattry, 0, (loops+1)*sizeof(double));
		SigmaJKtry = 0.0;
		SigmaJKspecheattry = 0.0;
		meanSpecheat = 0.0;

		numd = loops/d;
		for(j=0; j<numd; j++)
		{
			memset(EnergyJK, 0, (loops+1)*sizeof(double));
			memset(EnergyJK2, 0, (loops+1)*sizeof(double));
			//summing over omitted elements
			for(i=0; i<d; i++)
			{
				EnergyJK[j] += Data[d*j+i];
				EnergyJK2[j] += Data[d*j+i]*Data[d*j+i];
			}

			//subtracting omitted elements from total sum of energies
			EnergyJK[j] = (totEnergy - EnergyJK[j]) /(loops-d*1.0);
			//subtracting square of omitted elements from tot sum of energies for spec.heat
			EnergyJK2[j] = (totEnergy*totEnergy - EnergyJK2[j]) /(loops-d*1.0);
			//computing corrected square of sigma
			SigmaJKtry +=  (numd - 1.0)*(EnergyJK[j] - Mean)*(EnergyJK[j] - Mean) /numd;

			//computing specheat without omitted elements
			Specheattry[j] = (EnergyJK2[j] - EnergyJK[j]*EnergyJK[j])/numd;
			//computing mean specheat
			meanSpecheat += Specheattry[j]/numd;
		}
		for(j=0;j<numd;j++)
		{
			SigmaJKspecheattry += (numd - 1.0)*(Specheattry[j]-meanSpecheat)*(Specheattry[j]-meanSpecheat)/numd;
		}
		//Jackknife corrected standard deviation
		SigmaJKtry = sqrt(SigmaJKtry);
		if(SigmaJKtry > *SigmaJK)
		{
			//saving largest value of sigma corresponding to specific # of omitted elements
			*SigmaJK = SigmaJKtry;
			//Jackknife corrected standard deviation for specheat
			*SigmaJKspecheat = sqrt(SigmaJKspecheattry);
		}
	}

	return 0;
}

/*
 * computes integrated AutoCorrelation time
 */
double computeAutocorrelationtime(double *Data, double Mean, int loops, double Sigmast2)
{
	int i,d;
	double Autoc, Tau = 0.5;

	for(d=1; d<(loops/16); d++)
	{
		Autoc=0;
		for(i=0; i<(loops-d); i++)
		{
			Autoc += (Data[i]-Mean)*(Data[i+d]-Mean);
		}
		Autoc /= (Sigmast2*(loops*1.0-(d*1.0)));
		Tau += Autoc;
	}

	return Tau;
}
