#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "prng.h"
#include "montecarlo.h"
#include "hmc.h"
#include "EV.h"
#include "RepsSUn.h"
#include "model.h"
#include "MatrixMan.h"
#include "gammamatrices.h"
#include "RandomGens.h"

/*
 * code performs a Hybrid-Monte-Carlo simulation of Yang-Mills-Myers model for
 * 			-) 3 matrices: all files that include 3MM
 * 			-) 8 matrices: all files that include 8MM
 * 		and for modified Yang-Mills-Myers model for 8 matrices: also includes file S1.c, deltaS1.c
 * 		as a test code also includes Gaussian model: action.c
 * 	which model to a run the simulation for is specified in model.h
 */

// struct to measure length of simulation
struct tm tim;
// defines binwidth that is used for eigenvalue distributions
#define BINWIDTH	0.01
#define BINWIDTHDIRAC	0.05

#include <sys/types.h>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, MATRIX_SIZE=1, THERM=1, STEPS=1, STEP=1;
	double ALPHATILDE=4, EPS=1;
	double LOWERACC=70, UPPERACC=90, LAMBDA=100;
	double BINWIDTHCOMM=0.01;

	int i,j,k,t=0;

	int N;
	doublecomplex *pmat[NUMMAT], *pmatold[NUMMAT], *pmom[NUMMAT];
	double S=0, H=0, P=0, deltaH=0, Send=0, Hend=0, Pend=0, trYM1=0, trCS1=0, meanYM=0, meanCS=0, trYMold=0, trCSold=0;

	int acc=0;
	int acccheck=0;
	double mcrand=0;

	int ratio=0;
	int percentage=0;

	double expdH=0, meanH=0, meanS=0;

	double *SmallEV, *LargEV;
	int *numev, **pSavEVDist, *BinsEV, compEV=1;
	double **SavEVStore, *LargeEVStore, *SmallEVStore;
	int *numevStore;

	double ksi=0.1931833;
	int L;

	double C2=1, C3;
	double S0=0,S1=0;

	char pathout[1000], name[256], name2[256], name3[256], name4[256], name5[256];
	int conf=1, EVstore, numrep=0, *rep, diag=0;
	double mommult, prop1=1, prop2=1;

	doublecomplex *ptmp[NUMMAT], *ptmp2[NUMMAT];

	doublecomplex *pgammamult, **gammamat;
	int gammamultsize, numgammamult=0;
	double Smallgamma=0, Largegamma=0;
	int *SavGammaDist, BinsGamma=0, gammasize=0, compC=0;

	doublecomplex *pComm, *piComm;
	double SmallEVComm=1000000000000000, LargEVComm=-1000000000000000;
	int numevComm=0;
	int *SavEVCommDist, BinsEVComm=0, compComm=1;
	double *SavEVCommStore, LargeEVCommStore=0, SmallEVCommStore=0;
	int numevCommStore=0;

	doublecomplex *pComm2, *piComm2;
	double SmallEVComm2=1000000000000000, LargEVComm2=-1000000000000000;
	int numevComm2=0;
	int *SavEVCommDist2, BinsEVComm2=0, compComm2=0;
	double *SavEVCommStore2, LargeEVCommStore2=0, SmallEVCommStore2=0;
	int numevCommStore2=0;

	doublecomplex *pDirac;
	int numDirac=0, diracsize=0;
	double SmallDirac=0, LargeDirac=0;
	int *SavDiracDist, BinsDirac=0, compD=0;

	doublecomplex **ppLieGen, *pB;
	int BSize, NumEVB=0;
	double SmallEVB=0, LargeEVB=0;
	int *pSavEVBDist, BinsEVB=0, compB=0;
	double *SavEVBStore, LargeEVBStore=0, SmallEVBStore=0;
	int numevBStore=0;

	double trX2=0, meantrX2=0;
	int therm=0;

	int parsnum=2;
	char *parsname[2];
	double parsval[2];

	FILE *out;
	FILE *out2;
//	FILE *out3;
	FILE *out4;
	FILE *out5;
	char filename[256];
	char filename2[256];
//	char filename3[256];
	char filename4[256];
	char filename5[256];

	clock_t begin=0, end=0;
	int hours,min,sec;

	// loading parameters for simulation from file or command line
	if(argc>1)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &ALPHATILDE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, &LAMBDA, pathout, &conf,
				&EVstore, &mommult, &numrep, &rep, &diag, &prop1, &prop2, &compEV, &compComm, &compC, &compD, &compB);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &ALPHATILDE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, &LAMBDA, pathout, &conf, &mommult,
				&EVstore, &numrep, &rep, &diag, &prop1, &prop2, &compEV, &compComm, &compC, &compD, &compB);
	}
	sprintf(name, "%d", NUMMAT);
	strncat(name, "MMEV", 4);
	sprintf(name2, "%d", NUMMAT);
	strncat(name2, "MMCEV", 5);
	sprintf(name3, "%d", NUMMAT);
	strncat(name3, "MMCommEV", 8);
	sprintf(name4, "%d", NUMMAT);
	strncat(name4, "MMDiracEV", 9);
	sprintf(name5, "%d", NUMMAT);
	strncat(name5, "MMBEV", 5);

	sprintf(filename, "%s%dMMHMCLogfile-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, MATRIX_SIZE);
	sprintf(filename2, "%s%dMMHMCEnergy-%dx%d-a=%.4lf-l=%.2f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
//	sprintf(filename3, "%sHMCeps-%dx%d-a=%.2lf-l=%.2f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
	sprintf(filename4, "%s%dMMtrX2-%dx%d-a=%.2lf-l=%.4f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
	sprintf(filename5, "%s%dMMHMCYMCSS0S1-%dx%d-a=%.4lf-l=%.2f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
	out = fopen(filename, "a");
	out2 = fopen(filename2, "w");
//	out3 = fopen(filename3, "w");
	out4 = fopen(filename4, "w");
	out5 = fopen(filename5, "w");

	// start clock for time measurement
	begin=clock();
	time_t m_startTime;
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);

	// initializing seed for PRNG
	mt_goodseed();

	// printing parameters for simulation on screen and into logfile
	printf("Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	printf("%dMM code in use\n", NUMMAT);
	printf("loopnumber=%d\tsize=%d\tatilde=%lf\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\tlambda=%f start=%d\n",
			LOOPNUMBER, MATRIX_SIZE, ALPHATILDE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP, LAMBDA, conf);
	printf("pathout=%s\n", pathout);

	fprintf(out,"Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	fprintf(out,"%dMM HMC Simulation with %d Loops, %d integration steps, saving every %d step\n "
			"%dx%d matrices, Alphatilde = %lf, Thermalization = %d loops, starting with Epsilon = %lf, Range of Acceptance Rate is %lf to %lf\n"
			"%d integration loops, lambda = %lf\n",
			NUMMAT, LOOPNUMBER, STEPS, STEP, MATRIX_SIZE, MATRIX_SIZE, ALPHATILDE, THERM, EPS, UPPERACC, LOWERACC, STEPS, LAMBDA);
	fclose(out);

	N = MATRIX_SIZE*MATRIX_SIZE;
	gammasize = computeGammasize(NUMMAT);

	gammamultsize = MATRIX_SIZE*gammasize;
	diracsize = MATRIX_SIZE*MATRIX_SIZE*gammasize;
	BSize = LIEGROUP*MATRIX_SIZE;
	pgammamult = (doublecomplex*) malloc(gammamultsize*gammamultsize*sizeof(doublecomplex));
	gammamat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pDirac = (doublecomplex*) malloc(diracsize*diracsize*sizeof(doublecomplex));
	numev = (int*) malloc(sizeof(int)*NUMMAT);
	SmallEV = (double*) malloc(sizeof(double)*NUMMAT);
	LargEV = (double*) malloc(sizeof(double)*NUMMAT);
	pSavEVDist = (int**) malloc(sizeof(int*)*NUMMAT);
	BinsEV = (int*) malloc(sizeof(int)*NUMMAT);
	pComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	piComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	SavEVCommDist = (int*) malloc(sizeof(int*));
	pComm2 = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	piComm2 = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	SavEVCommDist2 = (int*) malloc(sizeof(int*));
	ppLieGen = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pB = (doublecomplex*) malloc(sizeof(doublecomplex)*BSize*BSize);
	SavEVStore = (double**) malloc(sizeof(double)*NUMMAT);
	LargeEVStore = (double*) malloc(sizeof(double)*NUMMAT);
	SmallEVStore = (double*) malloc(sizeof(double)*NUMMAT);
	numevStore = (int*) malloc(sizeof(int)*NUMMAT);
	SavEVCommStore = (double*) malloc(sizeof(double)*5000*MATRIX_SIZE);
	SavEVCommStore2 = (double*) malloc(sizeof(double)*5000*MATRIX_SIZE);
	SavEVBStore = (double*) malloc(sizeof(double)*5000*MATRIX_SIZE);
	for(i=0;i<parsnum;i++)
		parsname[i] = (char*) malloc(sizeof(char)*256);
	for(i=0;i<NUMMAT;i++)
	{
		pmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		ptmp[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		ptmp2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmatold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmom[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		gammamat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*gammasize*gammasize);
		ppLieGen[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*LIEGROUP*LIEGROUP);
		SavEVStore[i] = (double*) malloc(sizeof(double)*5000*MATRIX_SIZE);
	}

	memset(numev, 0, sizeof(int)*NUMMAT);
	memset(SmallEV, 0, sizeof(double)*NUMMAT);
	memset(LargEV, 0, sizeof(double)*NUMMAT);
	memset(BinsEV, 0, sizeof(int)*NUMMAT);
	memset(numevStore, 0, sizeof(int)*NUMMAT);
	memset(SmallEVStore, 0, sizeof(double)*NUMMAT);
	memset(LargeEVStore, 0, sizeof(double)*NUMMAT);
	memset(SavEVCommStore, 0, sizeof(double)*5000*MATRIX_SIZE);
	memset(SavEVCommStore2, 0, sizeof(double)*5000*MATRIX_SIZE);
	memset(SavEVBStore, 0, sizeof(double)*5000*MATRIX_SIZE);
	for(i=0;i<NUMMAT;i++)
	{
		memset(pmat[i], 0, sizeof(doublecomplex)*N);
		memset(ptmp[i], 0, sizeof(doublecomplex)*N);
		memset(ptmp2[i], 0, sizeof(doublecomplex)*N);
		memset(pmatold[i], 0, sizeof(doublecomplex)*N);
		memset(pmom[i], 0, sizeof(doublecomplex)*N);
		memset(gammamat[i], 0, sizeof(doublecomplex)*gammasize*gammasize);
		memset(ppLieGen[i], 0, sizeof(doublecomplex)*LIEGROUP*LIEGROUP);
		memset(SavEVStore[i], 0, sizeof(double)*5000*MATRIX_SIZE);
	}

	parsname[0] = "a";
	parsname[1] = "l";
	parsval[0] = ALPHATILDE;
	parsval[1] = LAMBDA;

	// obtain lie generators and gammamatrices for correct dimensionality (d=3 or d=8)
	RepsSUn(ppLieGen, LIEGROUP, LIEGROUP, NUMMAT, 1);
	gammamatrices(gammamat, NUMMAT);

	// define starting configuration
	startconf(pmat, ALPHATILDE, NUMMAT, MATRIX_SIZE, LIEGROUP, numrep, rep, conf, &L, diag, prop1, prop2);

	// compute Casimir operators; used for modified 8MM
	C2=(LIEGROUP-1.0)/(2.0*LIEGROUP) * L*(L+LIEGROUP);
	C3=/*C2**/((LIEGROUP*1.0-2.0)*(2.0*L+LIEGROUP*1.0))/(2.0*LIEGROUP);

	// compute value of action for current model
	S = action_function(pmat, &trYM1, &trCS1, ALPHATILDE, MATRIX_SIZE, NUMMAT, C2, C3, LAMBDA, &S0, &S1);
	trYMold=trYM1; trCSold=trCS1;
	printf("S=%f\tS0=%f\tS1=%f\n", S,S0,S1);

	// starting MC-loops
	for(k=1;k<LOOPNUMBER;k++)
	{
		//prints percentage on screen
		ratio = (k*1.0)/(LOOPNUMBER*1.0) *100.0;
		if(percentage<ratio)
		{
			percentage++;
			if(percentage%5==0)
				printf("%d", percentage);
			else if(percentage==100)
				printf("100\n");
			else
				printf("*");
			fflush(stdout);
		}

		//dynamical fit of interval for fluctuation
		if((k%100)==0)
		{
			if(acccheck < LOWERACC)
			{
				EPS = EPS*0.8;
			}
			else if(acccheck > UPPERACC)
			{
				EPS = EPS * 1.2;
			}
			acccheck=0;
		}

		// copying current configuration to pmatold if new conf not accepted
		for(i=0;i<NUMMAT;i++)
		{
			memcpy(pmatold[i], pmat[i], sizeof(doublecomplex)*N);
		}

		// generating gaussian momentum
		genmom(pmom, MATRIX_SIZE, mommult);

		// computing initial momentum and hamiltonian
		P = compmom(pmom, NUMMAT, MATRIX_SIZE, LAMBDA);
		H = S + P;
		trYM1=0;trCS1=0;

		// following hamiltonian equations of motion for number of steps = STEPS
		for(i=0;i<STEPS;i++)
		{
			//Leapfrog
//			addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);
//			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS, C2, C3, LAMBDA);
//			addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);

			//Omelyan
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, LAMBDA);
			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS/2.0, C2, C3, LAMBDA);
			addmat(pmat, pmom, EPS*(1-2*ksi), NUMMAT, MATRIX_SIZE, LAMBDA);
			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS/2.0, C2, C3, LAMBDA);
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, LAMBDA);
		}

		// computation of final value of action, momentum and hamiltonian
		Send = action_function(pmat, &trYM1, &trCS1, ALPHATILDE, MATRIX_SIZE, NUMMAT, C2, C3, LAMBDA, &S0, &S1);
		Pend = compmom(pmom, NUMMAT, MATRIX_SIZE, LAMBDA);
		Hend = Send + Pend;

		// comparing initial and final hamiltonian
		deltaH = H - Hend;

		//  ---------------- Metropolis step ---------------------------------
		if(deltaH>0)
		{
			H = Hend;
			S = Send;
			P = Pend;
			trYMold=trYM1; trCSold=trCS1;
			acc++;
			acccheck++;
		}
		else
		{
			mcrand = mt_ldrand();
			if(mcrand < exp(deltaH))
			{
				H = Hend;
				S = Send;
				P = Pend;
				trYMold=trYM1; trCSold=trCS1;
				acc++;
				acccheck++;
			}
			else
			{
				trYM1=trYMold; trCS1=trCSold;
				for(i=0;i<NUMMAT;i++)
				{
					memcpy(pmat[i], pmatold[i], sizeof(doublecomplex)*N);
				}
			}
		}
		fprintf(out2,"%f\n", S);

		// computation of observables every STEP steps if loopnumber is larger thermalization threshold
		if(k%STEP == 0 && k > THERM)
		{
			// save different terms of action in file
			fprintf(out5, "%f\t%f\t%f\t%f\n",trYM1, trCS1, S0, S1);
			fflush(out2);

			// compute code-checker exp(\delta H); should be close to 1
			expdH += exp(deltaH);

			// compute EV dist for matrices and save in pSavEVDist or store EV's in array SavEVStore
			if(compEV!=0)
			{
				if(EVstore==0){
					for(i=0;i<NUMMAT;i++)
						EVonthefly(&pSavEVDist[i], pmat[i], MATRIX_SIZE, BINWIDTH, &BinsEV[i], &LargEV[i], &SmallEV[i], &numev[i]);
				}
				else{
					for(i=0;i<NUMMAT;i++)
						SaveEV(SavEVStore[i], &LargeEVStore[i], &SmallEVStore[i], pmat[i], MATRIX_SIZE, &numevStore[i]);
				}
			}
			// compute EV dist for commutator i[X_1,X_2] or store EV's
			if(compComm!=0)
			{
				memset(pComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(piComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				Comm(pComm, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
				AddMat4(piComm, pComm, MATRIX_SIZE, 1.0);
				if(EVstore==0)
					EVonthefly(&SavEVCommDist, piComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm, &LargEVComm, &SmallEVComm, &numevComm);
				else
					SaveEV(SavEVCommStore, &LargeEVCommStore, &SmallEVCommStore, piComm, MATRIX_SIZE, &numevCommStore);
			}
			// compute EV dist for commutator i[X_1,X_8] or store EV's
			if(compComm2!=0){
				memset(pComm2, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(piComm2, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				Comm(pComm2, pmat[0], pmat[7], MATRIX_SIZE, 1.0);
				AddMat4(piComm2, pComm2, MATRIX_SIZE, 1.0);
				if(EVstore==0)
					EVonthefly(&SavEVCommDist, piComm2, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm2, &LargEVComm2, &SmallEVComm2, &numevComm2);
				else
					SaveEV(SavEVCommStore2, &LargeEVCommStore2, &SmallEVCommStore2, piComm2, MATRIX_SIZE, &numevCommStore2);
			}
			// compute EV distribution for matrix B = \sigma_{\mu} \otimes X_{\mu} where \mu=1,...,d
			if(compB!=0)
			{
				memset(pB, 0, sizeof(doublecomplex)*BSize*BSize);
				for(i=0;i<NUMMAT;i++)
				{
					Multilambda(pB, pmat[i], ppLieGen[i], MATRIX_SIZE, LIEGROUP);
				}
				if(EVstore==0)
					EVonthefly(&pSavEVBDist, pB, BSize, BINWIDTH, &BinsEVB, &LargeEVB, &SmallEVB, &NumEVB);
				else
					SaveEV(SavEVBStore, &LargeEVBStore, &SmallEVBStore, pB, BSize, &numevBStore);
			}
			// compute EV dist of matrix C = \gamma_{\mu} \otimes X_{\mu} where \mu = 1,...,d
			if(compC!=0)
			{
				memset(pgammamult, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
				for(i=0;i<NUMMAT;i++)
				{
					Multilambda(pgammamult, pmat[i], gammamat[i], MATRIX_SIZE, gammasize);
				}
				EVonthefly(&SavGammaDist, pgammamult, gammamultsize, BINWIDTH, &BinsGamma, &Largegamma, &Smallgamma, &numgammamult);
			}
			// compute EV dist of Dirac operator D = \gamma_{\mu} \otimes [X_{\mu}, . ] where \mu = 1,..,d
			if(compD!=0)
			{
				memset(pDirac, 0, sizeof(doublecomplex)*diracsize*diracsize);
				genD(pDirac, pmat, gammamat, MATRIX_SIZE, gammasize, NUMMAT);
				EVonthefly(&SavDiracDist, pDirac, diracsize, BINWIDTHDIRAC, &BinsDirac, &LargeDirac, &SmallDirac, &numDirac);
			}

			// compute \sum_{\mu} Tr(X_{\mu}^2)
			trX2=0;
			for(i=0;i<NUMMAT;i++)
			{
				diagMulti(&trX2, pmat[i], pmat[i], MATRIX_SIZE, 1.0);
			}
			meantrX2+=trX2;
			fprintf(out4, "%f\n",  trX2);
			fflush(out4);

			therm++;

			/* Save EV distributions every 5000th step */
			if((therm%5000)==0)
			{
				if(compEV!=0)
				{
					for(i=0;i<NUMMAT;i++)
					{
						FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
						if(EVstore==0)
							SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
						else
							StoreData(SavEVStore[i], &numevStore[i]);
						FileClose(EVstore);
					}
				}
				if(compComm!=0)
				{
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
					if(EVstore==0)
						SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
					else
						StoreData(SavEVCommStore, &numevCommStore);
					FileClose(EVstore);
				}
				if(compComm2!=0)
				{
					FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
					if(EVstore==0)
						SaveEVendofflight(SavEVCommDist2, BinsEVComm2, BINWIDTHCOMM, SmallEVComm2, numevComm2);
					else
						StoreData(SavEVCommStore2, &numevCommStore2);
					FileClose(EVstore);
				}
				if(compB!=0)
				{
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
					if(EVstore==0)
						SaveEVendofflight(pSavEVBDist, BinsEVB, BINWIDTH, SmallEVB, NumEVB);
					else
						StoreData(SavEVBStore, &numevBStore);
					FileClose(EVstore);
				}
				if(NUMMAT<14)
				{
					if(compC!=0)
					{
						FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name2, parsnum, parsname, parsval);
						SaveEVendofflight(SavGammaDist, BinsGamma, BINWIDTH, Smallgamma, numgammamult);
						FileClose(EVstore);
					}
					if(compD!=0)
					{
						FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name4, parsnum, parsname, parsval);
						SaveEVendofflight(SavDiracDist, BinsDirac, BINWIDTHDIRAC, SmallDirac, numDirac);
						FileClose(EVstore);
					}
				}
			}
			meanH += H; meanS += S;
			meanYM += trYM1; meanCS += trCS1;
		}
	}
	printf("\n");
	// save final EV distribution results
	if(compEV!=0)
	{
		for(i=0;i<NUMMAT;i++)
		{
			FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
			if(EVstore==0)
				SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
			else
				StoreData(SavEVStore[i], &numevStore[i]);
			FileClose(EVstore);
		}
	}
	if(compComm!=0)
	{
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
		if(EVstore==0)
			SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
		else
			StoreData(SavEVCommStore, &numevCommStore);
		FileClose(EVstore);
	}
	if(compComm2!=0)
	{
		FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
		if(EVstore==0)
			SaveEVendofflight(SavEVCommDist2, BinsEVComm2, BINWIDTHCOMM, SmallEVComm2, numevComm2);
		else
			StoreData(SavEVCommStore2, &numevCommStore2);
		FileClose(EVstore);
	}
	if(compB!=0)
	{
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
		if(EVstore==0)
			SaveEVendofflight(pSavEVBDist, BinsEVB, BINWIDTH, SmallEVB, NumEVB);
		else
			StoreData(SavEVBStore, &numevBStore);
		FileClose(EVstore);
	}
	if(NUMMAT<14)
	{
		if(compC!=0)
		{
			FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name2, parsnum, parsname, parsval);
			SaveEVendofflight(SavGammaDist, BinsGamma, BINWIDTH, Smallgamma, numgammamult);
			FileClose(EVstore);
		}
		if(compD!=0)
		{
			FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name4, parsnum, parsname, parsval);
			SaveEVendofflight(SavDiracDist, BinsDirac, BINWIDTHDIRAC, SmallDirac, numDirac);
			FileClose(EVstore);
		}
	}

	// print results on screen and into logfile
	end=clock();
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	hours= (int) (((1.*(end - begin)) / CLOCKS_PER_SEC)/60)/60;
	min= (int)((1.*(end - begin)) / CLOCKS_PER_SEC)/60 - (int)((((1.*(end - begin)) / CLOCKS_PER_SEC)/60)/60)*60;
	sec= (int)((1.*(end - begin)) / CLOCKS_PER_SEC) - (int)(((1.*(end - begin)) / CLOCKS_PER_SEC)/60)/60 * 3600 - (int)(((1.*(end - begin)) / CLOCKS_PER_SEC)/60)*60;
	if(LIEGROUP==2)
		printf("SU(2) ground state = %f\n", -N*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12);
	if((LIEGROUP-2)!=0)
	{
		printf("SU(3) ground state = %f\tSU(2) ground state = %f\n",
				-MATRIX_SIZE*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12,
				-MATRIX_SIZE*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*(N-1)/24);
	}
	printf("acc=%d\tacc rate = %.2lf\n", acc, 100.0*acc/(k*1.0));
	printf("<expdH>=%f\t<H>=%f\t<S>=%f\n", expdH/(therm*1.0), meanH/(therm*1.0), meanS/(therm*1.0));
	printf("ID = %f\t<trX2/N>=%f\n", (4*meanYM + 3*meanCS)/(therm*1.0*(N-1.0)), trX2/(1.*therm*MATRIX_SIZE));
	printf("Finished at %d:%d:%d\t%d/%d\nTime elapsed: %d h:%d min:%d s \n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1,hours, min,sec);

	out = fopen(filename, "a");
	fprintf(out,"<expdH> = %f\tAcceptance Rate = %f\n<S> = %f\n",
			expdH/(t*1.0), 100.0*acc/(k*1.0), meanS/(t*1.0));
	if(LIEGROUP==2)
		fprintf(out, "SU(2) ground state = %f\n", -N*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12);
	if((LIEGROUP-2)!=0)
	{
		fprintf(out, "SU(3) ground state = %f\tSU(2) ground state = %f\n",
				-MATRIX_SIZE*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12,
				-MATRIX_SIZE*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*(N-1)/24);
	}
	fprintf(out, "ID = %f\t<trX2/N>=%f\n", (4*meanYM + 3*meanCS)/(therm*1.0*(N-1.0)), trX2/(1.*therm*MATRIX_SIZE));
	fprintf(out,"Finished at %d:%d:%d\t%d/%d\nTime elapsed: %d h:%d min:%d s \n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1,hours, min,sec);

	fclose(out);
	fclose(out2);
//	fclose(out3);
	fclose(out4);
	fclose(out5);

	for(i=0;i<NUMMAT;i++)
	{
		free(pmat[i]);
		free(pmatold[i]);
		free(pmom[i]);
		free(gammamat[i]);
		free(ppLieGen[i]);
	}
	free(ppLieGen);
	free(pB);
	free(pComm);
	free(piComm);
	free(pComm2);
	free(piComm2);
	free(pgammamult);

	return 0;
}
