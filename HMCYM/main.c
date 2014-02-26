#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"
#include "EV.h"
#include "RepsSUn.h"
#include "MatrixMan.h"
#include "gammamatrices.h"
#include "RandomGens.h"

/*
 * Hybrid Monte-Carlo of pure Yang-Mills model with D matrices
 */

struct tm tim;
FILE *outTraces;

#include <sys/types.h>
#include <sys/stat.h>

#define BINWIDTHDIRAC	0.05

int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, MATRIX_SIZE=1, THERM=1, STEPS=1, STEP=1, NUMMAT=3;
	double EPS=1, LOWERACC=70, UPPERACC=90, BINWIDTH=0.01, BINWIDTHCOMM=0.01;

	int i,k;

	int N;
	doublecomplex **pmat, **pmatold, **pmom;
	double S=0, H=0, P=0, deltaH=0, Send=0, Hend=0, Pend=0;

	int acc=0;
	int acccheck=0;
	double mcrand=0;

	int ratio=0;											//for percentage on screen
	int percentage=0;

	double expdH=0, meanH=0, meanS=0;

	double *SmallEV, *LargEV;
	int *numev, **pSavEVDist, *BinsEV, compEV=1;
	double ksi=0.1931833;

	char pathout[1000], name[256], name2[256], name3[256], name4[256], name5[256];
	int parsnum=0;
	char *parsname[1];
	double parsval[1];
	int EVstore;

	int therm=0;

	doublecomplex *pgammamult;
	int gammamultsize, numgammamult=0;
	double Smallgamma=0, Largegamma=0;
	int *SavGammaDist, BinsGamma=0;
	int gammasize=0, compC=0;
	doublecomplex **gammamat;

	doublecomplex *pComm, *piComm;
	double SmallEVComm=1000000000000000, LargEVComm=-1000000000000000;
	int numevComm=0, compComm=1;
	int *SavEVCommDist, BinsEVComm=0;

	doublecomplex *pDirac;
	int numDirac=0, diracsize=0;
	double SmallDirac=0, LargeDirac=0;
	int *SavDiracDist, BinsDirac=0, compD=0;

	doublecomplex **tmpmat;
	double trXi2=0, meantrXi2=0;
	int compTraces=0, compTraceX2=0;

	int compGdd=0;
	doublecomplex *matGdd;
	double largeEVGdd=-100000, smallEVGdd=10000;
	int *saveEVGdd, numEVGdd=0, binsGdd=0;

	double mass=0;
	int start=0;

	FILE *out;
	FILE *out2;
	FILE *out3;
	FILE *out4;
//	FILE *out5;
//	FILE *out6;
	char filename[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
//	char filename5[256];
//	char filename6[256];
	char filenameTr[256];

	// loading parameters for simulation from file or command line
	if(argc>1)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &NUMMAT,
				&mass, &start, &compC, &compD, &compEV, &compComm, &compTraces, &compTraceX2, &compGdd);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &NUMMAT, &mass,
				&start, &compC, &compD, &compEV, &compComm, &compTraces, &compTraceX2, &compGdd);
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
	strncat(name5, "MMGddEV", 7);

	sprintf(filename, "%sYM%dMMHMCLogfile-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, MATRIX_SIZE);
	if(mass==0)
	{
		BINWIDTH = 0.01;BINWIDTHCOMM=0.01;
		sprintf(filename2, "%sYM%dMMHMCEnergy-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
		sprintf(filename3, "%sHMCeps-%dx%d.txt", pathout, MATRIX_SIZE, LOOPNUMBER);
		sprintf(filename4, "%s%dMMtrXYZ-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
//		sprintf(filename5, "%sYM%dMMHMCEnergyS0+S1-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
//		sprintf(filename6, "%sYM%dMMlambdaMatrixEV-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
	}
	if(mass!=0)
	{
		BINWIDTH = 0.002;BINWIDTHCOMM=0.0002;
		sprintf(filename2, "%sYM%dMMHMCEnergy-%dx%d-m=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, mass);
		sprintf(filename3, "%sHMCeps-%dx%d-m=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, mass);
		sprintf(filename4, "%s%dMMtrXYZ-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
//		sprintf(filename5, "%sYM%dMMHMCEnergyS0+S1-%dx%d-m=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, mass);
//		sprintf(filename6, "%sYM%dMMlambdaMatrixEV-%dx%d-m=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, mass);
	}
	if(compTraces!=0)
	{
		sprintf(filenameTr, "%s%dMMtraces68-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
		outTraces = fopen(filenameTr, "w");
	}
	out = fopen(filename, "a");
	out2 = fopen(filename2, "w");
	out3 = fopen(filename3, "w");
	out4 = fopen(filename4, "w");
//	out5 = fopen(filename5, "w");
//	out6 = fopen(filename6, "w");

	// initializing seed for PRNG
	mt_goodseed();

	// printing parameters for simulation on screen and into logfile
	time_t m_startTime;
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	printf("YM with %d matrices\n", NUMMAT);
	printf("loopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\n",
			LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP);
	printf("pathout=%s\n", pathout);

	fprintf(out,"Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	fprintf(out,"YM HMC Simulation with %d matrices for %d loops, %d integration steps, saving every %d step\n "
			"%dx%d matrices, Thermalization = %d loops, starting with Epsilon = %lf, Range of Acceptance Rate is %lf to %lf\n"
			"%d integration loops\n",
			NUMMAT, LOOPNUMBER, STEPS, STEP, MATRIX_SIZE, MATRIX_SIZE, THERM, EPS, UPPERACC, LOWERACC, STEPS);
	fclose(out);

	N = MATRIX_SIZE*MATRIX_SIZE;

	for(i=0;i<parsnum;i++)
		parsname[i] = (char*) malloc(sizeof(char)*256);
	pmat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pmatold = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pmom = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	numev = (int*) malloc(sizeof(int)*NUMMAT);
	SmallEV = (double*) malloc(sizeof(double)*NUMMAT);
	LargEV = (double*) malloc(sizeof(double)*NUMMAT);
	pSavEVDist = (int**) malloc(sizeof(int*)*NUMMAT);
	BinsEV = (int*) malloc(sizeof(int)*NUMMAT);
	pComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	piComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	SavEVCommDist = (int*) malloc(sizeof(int*));
	tmpmat = (doublecomplex**) malloc(sizeof(doublecomplex)*3);
	matGdd = (doublecomplex*) malloc(sizeof(doublecomplex)*NUMMAT*NUMMAT);
	for(i=0;i<NUMMAT;i++)
	{
		pmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmatold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmom[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	}
	for(i=0;i<3;i++)
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);

	memset(numev, 0, sizeof(int)*NUMMAT);
	memset(SmallEV, 0, sizeof(double)*NUMMAT);
	memset(LargEV, 0, sizeof(double)*NUMMAT);
	memset(BinsEV, 0, sizeof(int)*NUMMAT);
	for(i=0;i<NUMMAT;i++)
	{
		memset(pmat[i], 0, sizeof(doublecomplex)*N);
		memset(pmatold[i], 0, sizeof(doublecomplex)*N);
		memset(pmom[i], 0, sizeof(doublecomplex)*N);
	}
	for(i=0;i<3;i++)
		memset(tmpmat[i], 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);

	// execute this code only if number of matrices NUMMAT is smaller 14; otherwise computation takes too long
	if(NUMMAT<14)
	{
		gammasize = computeGammasize(NUMMAT);
		gammamultsize = MATRIX_SIZE*gammasize;
		diracsize = N*gammasize;
		pgammamult = (doublecomplex*) malloc(gammamultsize*gammamultsize*sizeof(doublecomplex));
		gammamat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
		if(compD!=0)
		{
			pDirac = (doublecomplex*) malloc(diracsize*diracsize*sizeof(doublecomplex));
		}
		for(i=0;i<NUMMAT;i++)
		{
			gammamat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*gammasize*gammasize);
			memset(gammamat[i], 0, sizeof(doublecomplex)*gammasize*gammasize);
		}

		// generating \gamma-matrices corresponding to dimensionality of model; only necessary for computation of matric C and D
		gammamatrices(gammamat, NUMMAT);
		printf("gammasize=%d\tgammamultsize=%d\n", gammasize, gammamultsize);
	}

	// if start==1 start from random configuration; otherwise matrices are initialized to zero
	if(start==1)
		gen_randphicplx(pmat, NUMMAT, MATRIX_SIZE);

	// computation of value of initial action
	S = actionYM(pmat, MATRIX_SIZE, NUMMAT, mass);

	printf("S=%f\n", S);
//	return 0;

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
		fprintf(out3, "%f\n", EPS);

		// copying current configuration to pmatold if new conf not accepted
		for(i=0;i<NUMMAT;i++)
		{
			memcpy(pmatold[i], pmat[i], sizeof(doublecomplex)*N);
		}

		// generating gaussian momentum
		gen_gaussmomcplx(pmom, NUMMAT, MATRIX_SIZE, 1.0);

		// computing initial momentum
		P = mom(pmom, NUMMAT, MATRIX_SIZE);
		H = S + P;

		// following hamiltonian equations of motion for number of steps = STEPS
		for(i=0;i<STEPS;i++)
		{
			//Leapfrog
//			addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);
//			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS, C2, C3, LAMBDA);
//			addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);

			//Omelyan
			addmattrace(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE);
			deltaYM(pmat, pmom, NUMMAT, MATRIX_SIZE, EPS/2.0, mass);
			addmattrace(pmat, pmom, EPS*(1-2*ksi), NUMMAT, MATRIX_SIZE);
			deltaYM(pmat, pmom, NUMMAT, MATRIX_SIZE, EPS/2.0, mass);
			addmattrace(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE);
		}
		// computation of final value of action, momentum and hamiltonian
		Send = actionYM(pmat, MATRIX_SIZE, NUMMAT, mass);
		Pend = mom(pmom, NUMMAT, MATRIX_SIZE);
		Hend = Send + Pend;

		// comparing initial and final hamiltonian
		deltaH = H - Hend;

		//  ---------------- Metropolis step ---------------------------------
		if(deltaH>0)
		{
			H = Hend;
			S = Send;
			P = Pend;
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
				acc++;
				acccheck++;
			}
			else
			{
				for(i=0;i<NUMMAT;i++)
				{
					memcpy(pmat[i], pmatold[i], sizeof(doublecomplex)*N);
				}
			}
		}

		fprintf(out2,"%f\n", S);
		fflush(out2);

		// computation of observables every STEP steps if loopnumber is larger thermalization threshold
		if(k%STEP == 0 && k > THERM)
		{
			// compute code-checker exp(\delta H); should be close to 1
			expdH += exp(deltaH);

			// compute EV dist for matrices
			if(compEV!=0)
			{
				for(i=0;i<NUMMAT;i++)
					EVonthefly(&pSavEVDist[i], pmat[i], MATRIX_SIZE, BINWIDTH, &BinsEV[i], &LargEV[i], &SmallEV[i], &numev[i]);
			}
			// compute EV dist for commutator i[X_1,X_2]
			if(compComm!=0)
			{
				memset(pComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(piComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				Comm(pComm, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
				AddMat4(piComm, pComm, MATRIX_SIZE, 1.0);
				EVonthefly(&SavEVCommDist, piComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm, &LargEVComm, &SmallEVComm, &numevComm);
			}

			// if number of matrices is smaller 14 check for computation of following observables
			if(NUMMAT<14)
			{
				// compute EV dist of matrix C = \gamma_1 \otimes X_1 + \gamma_2 \otimes X_2
				if(compC!=0)
				{
					memset(pgammamult, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
					for(i=0;i<2;i++)
					{
						Multilambda(pgammamult, pmat[i], gammamat[i], MATRIX_SIZE, gammasize);
					}
					EVonthefly(&SavGammaDist, pgammamult, gammamultsize, BINWIDTH, &BinsGamma, &Largegamma, &Smallgamma, &numgammamult);
				}
				// compute EV dist of Dirac operator D = \gamma_{\mu} \otimes [X_{\mu}, . ]
				if(compD!=0)
				{
					memset(pDirac, 0, sizeof(doublecomplex)*diracsize*diracsize);
					genD(pDirac, pmat, gammamat, MATRIX_SIZE, gammasize, NUMMAT);
					EVonthefly(&SavDiracDist, pDirac, diracsize, BINWIDTHDIRAC, &BinsDirac, &LargeDirac, &SmallDirac, &numDirac);
//					printf("numDirac=%d\tSmallDirac=%f\tLargeDirac=%f\n", numDirac, SmallDirac, LargeDirac);
				}
			}
			// compute \sum_{\mu} Tr(X_{\mu}^2)
			if(compTraceX2!=0){
				trXi2=0;
				for(i=0;i<NUMMAT;i++)																					//Radius Tr(X_aX^a)
				{
					diagMulti(&trXi2, pmat[i], pmat[i], MATRIX_SIZE, 1.0);
				}
			}
			// compute further correlators
			if(compTraces!=0)
				Traces(pmat, MATRIX_SIZE, NUMMAT, LOOPNUMBER, pathout);
			// compute EV dist of stress-energy tensor a la Nishimura; Gdd = Tr(X_{\mu}X_{\nu}) / N
			if(compGdd!=0){
				computeGdd(matGdd, pmat, NUMMAT, MATRIX_SIZE);
				EVonthefly(&saveEVGdd, matGdd, NUMMAT, BINWIDTH, &binsGdd, &largeEVGdd, &smallEVGdd, &numEVGdd);
			}

			therm++;

			/* Save EV distributions every 5000th step */
			if((therm%5000)==0)
			{
				if(compEV!=0)
				{
					for(i=0;i<NUMMAT;i++)
					{
						FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
						SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
						FileClose(EVstore);
					}
				}
				if(compComm!=0)
				{
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
					SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
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
				if(compGdd!=0){
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
					SaveEVendofflight(saveEVGdd, binsGdd, BINWIDTH, smallEVGdd, numEVGdd);
					FileClose(EVstore);
				}
			}
			meanH += H; meanS += S;
		}
	}
	printf("\n");
	// save final EV distribution results
	if(compEV!=0)
	{
		for(i=0;i<NUMMAT;i++)
		{
			FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
			SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
			FileClose(EVstore);
		}
	}
	if(compComm!=0)
	{
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
		SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
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
	if(compGdd!=0){
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
		SaveEVendofflight(saveEVGdd, binsGdd, BINWIDTH, smallEVGdd, numEVGdd);
		FileClose(EVstore);
	}

	// print results on screen and into logfile
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("YM with %d matrices\n", NUMMAT);
	printf("loopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\n",
			LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP);
	printf("pathout=%s\n", pathout);
	printf("acc=%d\tacc rate = %.2lf\n", acc, 100.0*acc/(k*1.0));
	printf("<expdH>=%f\t<trXi2/N>=%f\t<H>=%f\t<S>=%f\n", expdH/(therm*1.0), meantrXi2/(1.*therm*MATRIX_SIZE), meanH/(therm*1.0), meanS/(therm*1.0));
	printf("ID = %f\n", (4*meanS)/(therm*1.0*(N-1.0)));

	printf("Finished at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	out = fopen(filename, "a");
	fprintf(out,"<expdH> = %f\tAcceptance Rate = %f\n<S> = %f\t<trXi2/N>=%f\n",
			expdH/(therm*1.0), 100.0*acc/(k*1.0), meanS/(therm*1.0), meantrXi2/(1.*therm*MATRIX_SIZE));
	fprintf(out, "ID=%f\n", (4*meanS)/(therm*1.0*(N-1.0)));
	fprintf(out,"Finished at %d:%d:%d\t%d/%d\n\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	fclose(out);
	fclose(out2);
	fclose(out3);
	fclose(out4);
//	fclose(out5);
//	fclose(out6);

	for(i=0;i<NUMMAT;i++)
	{
		free(pmat[i]);
		free(pmatold[i]);
		free(pmom[i]);
	}
	free(pComm);
	free(piComm);
	if(NUMMAT<14)
	{
		for(i=0;i<NUMMAT;i++)
			free(gammamat[i]);
		free(pgammamult);
	}
	return 0;
}
