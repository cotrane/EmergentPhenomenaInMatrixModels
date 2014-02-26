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
#include "MatrixMan.h"
#include "gammamatrices.h"
#include "RandomGens.h"

struct tm tim;

#include <sys/types.h>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, MATRIX_SIZE=1, THERM=1, STEPS=1, STEP=1, NUMMAT=2, LIEGROUP=2;
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
	int *numev, **pSavEVDist, *BinsEV;
	double ksi=0.1931833;

	char pathout[1000], name[256], name2[256], name3[256], name5[256], name4[256];
	int parsnum=0;
	char *parsname[1];
	double parsval[1];
	int EVstore;

	double trX2=0, trX2mean=0;
	int therm=0;

	doublecomplex *pComm, *piComm;
	double SmallEVComm=1000000000000000, LargEVComm=-1000000000000000;
	int numevComm=0;
	int *SavEVCommDist, BinsEVComm=0;

	doublecomplex *pgammamult;
	int gammamultsize, numgammamult=0;
	double Smallgamma=0, Largegamma=0;
	int *SavGammaDist, BinsGamma=0;
	int gammasize=0, compC=0;
	doublecomplex **gammamat;

	int *pSavDistnHr, *pSavDistnHi, BinsEVnH[2];
	int numEVnH=0;
	doublecomplex LargEVnH, SmallEVnH, *pphi;
	int *pSavDistMod, BinsEVMod=0, numEVMod=0;
	double LargeEVMod=0, SmallEVMod=0;

	double mass=0;

	FILE *out;
	FILE *out2;
	FILE *out3;
	char filename[256];
	char filename2[256];
	char filename3[256];

	if(argc>1)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &NUMMAT,
				&LIEGROUP, &mass, &compC);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &NUMMAT, &LIEGROUP,
				&mass, &compC);
	}
	sprintf(name, "%d", NUMMAT);
	strncat(name, "MMEV", 4);
	sprintf(name2, "%d", NUMMAT);
	strncat(name2, "MMCEV", 5);
	sprintf(name3, "%d", NUMMAT);
	strncat(name3, "MMCommEV", 8);
	sprintf(name4, "%d", NUMMAT);
	strncat(name4, "MMmodEV", 7);
	sprintf(name5, "%d", NUMMAT);
	strncat(name5, "MMphiEV", 7);

	sprintf(filename, "%s%dMMHMCLogfile-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, MATRIX_SIZE);
	sprintf(filename2, "%s%dMMHMCEnergy-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
	sprintf(filename3, "%sHMCeps-%dx%d.txt", pathout, MATRIX_SIZE, LOOPNUMBER);

	out = fopen(filename, "a");
	out2 = fopen(filename2, "w");
	out3 = fopen(filename3, "w");

	time_t m_startTime;
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	printf("Free Gaussian with %d matrices\n", NUMMAT);
	printf("loopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\n",
			LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP);
	printf("pathout=%s\n", pathout);

	fprintf(out,"Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	fprintf(out,"Free Gaussian HMC Simulation with %d matrices for %d loops, %d integration steps, saving every %d step\n "
			"%dx%d matrices, Thermalization = %d loops, starting with Epsilon = %lf, Range of Acceptance Rate is %lf to %lf\n"
			"%d integration loops\n",
			NUMMAT, LOOPNUMBER, STEPS, STEP, MATRIX_SIZE, MATRIX_SIZE, THERM, EPS, UPPERACC, LOWERACC, STEPS);
	fclose(out);

	N = MATRIX_SIZE*MATRIX_SIZE;
	BINWIDTH = 0.01;BINWIDTHCOMM=0.01;

//	LargEVnH.r=0;LargEVnH.i=0;SmallEVnH.r=0;SmallEVnH.i=0;
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
	SavEVCommDist = (int*) malloc(sizeof(int));
	pphi = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	for(i=0;i<NUMMAT;i++)
	{
		pmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmatold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmom[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	}

	memset(numev, 0, sizeof(int)*NUMMAT);
	memset(SmallEV, 0, sizeof(double)*NUMMAT);
	memset(LargEV, 0, sizeof(double)*NUMMAT);
	memset(BinsEV, 0, sizeof(int)*NUMMAT);
	memset(pphi, 0, sizeof(doublecomplex)*N);
	memset(BinsEVnH, 0, sizeof(int)*2);
	for(i=0;i<NUMMAT;i++)
	{
		memset(pmat[i], 0, sizeof(doublecomplex)*N);
		memset(pmatold[i], 0, sizeof(doublecomplex)*N);
		memset(pmom[i], 0, sizeof(doublecomplex)*N);
	}
	if(NUMMAT<14 && NUMMAT>1 && compC!=0)
	{
		gammasize = computeGammasize(NUMMAT);
		gammamultsize = MATRIX_SIZE*gammasize;

		pgammamult = (doublecomplex*) malloc(gammamultsize*gammamultsize*sizeof(doublecomplex));
		gammamat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
		for(i=0;i<NUMMAT;i++)
		{
			gammamat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*gammasize*gammasize);
			memset(gammamat[i], 0, sizeof(doublecomplex)*gammasize*gammasize);
		}
		gammamatrices(gammamat, NUMMAT);
		printf("gammasize=%d\tgammamultsize=%d\n", gammasize, gammamultsize);
	}

	// generate random hermitian matrix
	gen_randphicplxtrace(pmat, NUMMAT, MATRIX_SIZE);

	// compute gaussian action
	S = action(pmat, NUMMAT, MATRIX_SIZE, mass);
	printf("S=%f\n", S);
//	return 0;

	for(k=1;k<LOOPNUMBER;k++)
	{
		// prints percentage on screen
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

		// dynamical fit of interval for fluctuation
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

		// copy current pmat into pmatold
		for(i=0;i<NUMMAT;i++)
		{
			memcpy(pmatold[i], pmat[i], sizeof(doublecomplex)*N);
		}

		// generate gaussian hermitian matrix pmom
		gen_gaussmomcplxtrace(pmom, NUMMAT, MATRIX_SIZE, 1.0);

//		printmat("pmom:", pmom[0], MATRIX_SIZE);

		// compute momentum and hamiltonian energy
		P = mom(pmom, NUMMAT, MATRIX_SIZE);
		H = S + P;

//		if(k==1)
//			fprintf(out5, "%f\t%f\t%f\n",H, S, P);

		// integrate over number of steps=STEPS
		for(i=0;i<STEPS;i++)
		{
			//Omelyan
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE);
			addmom(pmom, pmat, EPS/2.0, NUMMAT, MATRIX_SIZE, mass);
			addmat(pmat, pmom, EPS*(1-2*ksi), NUMMAT, MATRIX_SIZE);
			addmom(pmom, pmat, EPS/2.0, NUMMAT, MATRIX_SIZE, mass);
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE);
		}
		// compute final value of action, momentum and hamiltonian
		Send = action(pmat, NUMMAT, MATRIX_SIZE, mass);
		Pend = mom(pmom, NUMMAT, MATRIX_SIZE);
		Hend = Send + Pend;

//		fprintf(out5, "%f\t%f\t%f\n",Hend, Send, Pend);

		// compute difference between initial and final H
		deltaH = H - Hend;

		// accept if Hend < H
		if(deltaH>0)
		{
			H = Hend;
			S = Send;
			P = Pend;
			acc++;
			acccheck++;
		}
		// compare with random number if  Hend > H
		else
		{
			mcrand = prng_get_double();
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

		// compute Observables
		if(k%STEP == 0 && k > THERM)
		{
			// compute \sum e^{\Delta H}; should be \sim 1 when averaged over all mc-steps (averaging is done in the end)
			expdH += exp(deltaH);

			// compute Eigenvalues
			for(i=0;i<NUMMAT;i++)
				EVonthefly(&pSavEVDist[i], pmat[i], MATRIX_SIZE, BINWIDTH, &BinsEV[i], &LargEV[i], &SmallEV[i], &numev[i]);

			// compute <\sum_{\mu} Tr(pmat_{\mu}^2)>/MATRIX_SIZE; averaging is done in the end below
			trX2=0;
			for(i=0;i<NUMMAT;i++)
			{
				diagMulti(&trX2, pmat[i], pmat[i], MATRIX_SIZE, 1.0/MATRIX_SIZE);
			}
			trX2mean+=trX2;

			if(NUMMAT>1)
			{
				// compute matrix pphi = pmat[0] + i*pmat[1] and EV distribution of real and imaginary part as well as of modulus
				memset(pphi, 0, sizeof(doublecomplex)*N);
				AddMat9(pphi, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
				EVtoDist_nonHermitian(&pSavDistnHr, &pSavDistnHi, pphi, MATRIX_SIZE, BINWIDTH, BinsEVnH, &LargEVnH,
						&SmallEVnH, &numEVnH);
				EVtoDist_Modulus_nH(&pSavDistMod, pphi, MATRIX_SIZE, BINWIDTH, &BinsEVMod, &LargeEVMod, &SmallEVMod,
						&numEVMod);

				// compute i[pmatXYZ[0],pmatXYZ[1]] and its EV distribution
				memset(pComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(piComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				Comm(pComm, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
				AddMat4(piComm, pComm, MATRIX_SIZE, 1.0);
				EVonthefly(&SavEVCommDist, piComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm, &LargEVComm, &SmallEVComm, &numevComm);
			}

			// compute matrix C = \gamma_{\mu} \otimes pmatXYZ_{\mu}; \mu=1,2; and its EV distribution
			if(NUMMAT<14 && NUMMAT>1 && compC!=0)
			{
				memset(pgammamult, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
				for(i=0;i<2;i++)
				{
					Multilambda(pgammamult, pmat[i], gammamat[i], MATRIX_SIZE, gammasize);
				}
				EVonthefly(&SavGammaDist, pgammamult, gammamultsize, BINWIDTH, &BinsGamma, &Largegamma, &Smallgamma, &numgammamult);
			}

			meanH += H; meanS += S;
			therm++;

			/* Save EV distributions every 5000th step */
			if((therm%5000)==0)
			{
				for(i=0;i<NUMMAT;i++){
					FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
					SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
					FileClose(EVstore);
				}
				if(NUMMAT>1){
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
					SaveEVendofflight(pSavDistnHr, BinsEVnH[0], BINWIDTH, SmallEVnH.r, numEVnH);
					FileClose(EVstore);

					FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
					SaveEVendofflight(pSavDistnHi, BinsEVnH[1], BINWIDTH, SmallEVnH.i, numEVnH);
					FileClose(EVstore);

					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name4, parsnum, parsname, parsval);
					SaveEVendofflight(pSavDistMod, BinsEVMod, BINWIDTH, SmallEVMod, numEVMod);
					FileClose(EVstore);

					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
					SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
					FileClose(EVstore);
				}
				if(NUMMAT<14 && NUMMAT>1 && compC!=0){
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name2, parsnum, parsname, parsval);
					SaveEVendofflight(SavGammaDist, BinsGamma, BINWIDTH, Smallgamma, numgammamult);
					FileClose(EVstore);
				}
			}
		}
	}
	printf("\n");

	// save observables in the end
	for(i=0;i<NUMMAT;i++){
		FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
		SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
		FileClose(EVstore);
	}
	if(NUMMAT>1){
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
		SaveEVendofflight(pSavDistnHr, BinsEVnH[0], BINWIDTH, SmallEVnH.r, numEVnH);
		FileClose(EVstore);

		FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
		SaveEVendofflight(pSavDistnHi, BinsEVnH[1], BINWIDTH, SmallEVnH.i, numEVnH);
		FileClose(EVstore);

		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name4, parsnum, parsname, parsval);
		SaveEVendofflight(pSavDistMod, BinsEVMod, BINWIDTH, SmallEVMod, numEVMod);
		FileClose(EVstore);

		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
		SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
		FileClose(EVstore);
	}
	if(NUMMAT<14 && NUMMAT>1 && compC!=0){
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name2, parsnum, parsname, parsval);
		SaveEVendofflight(SavGammaDist, BinsGamma, BINWIDTH, Smallgamma, numgammamult);
		FileClose(EVstore);
	}

	// print some info to screen
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("YM with %d matrices\n", NUMMAT);
	printf("loopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\n",
			LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP);
	printf("pathout=%s\n", pathout);
	printf("acc=%d\tacc rate = %.2lf\n", acc, 100.0*acc/(k*1.0));
	printf("<expdH>=%f\t<H>=%f\t<S>=%f\t<trX2/N>=%f\n", expdH/(therm*1.0), meanH/(therm*1.0), meanS/(therm*1.0), trX2mean/(therm*1.));
	printf("ID = %f\n", (4*meanS)/(therm*1.0*(N-1.0)));
	printf("Finished at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	out = fopen(filename, "a");
	fprintf(out,"<expdH> = %f\tAcceptance Rate = %f\n<S> = %f\t <trX2>/N=%f\n",
			expdH/(therm*1.0), 100.0*acc/(k*1.0), meanS/(therm*1.0), trX2mean/(therm*1.*MATRIX_SIZE));
	fprintf(out, "ID=%f\n", (4*meanS)/(therm*1.0*(N-1.0)));
	fprintf(out,"Finished at %d:%d:%d\t%d/%d\n\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);


	fclose(out);
	fclose(out2);
	fclose(out3);

	for(i=0;i<NUMMAT;i++)
	{
		free(pmat[i]);
		free(pmatold[i]);
		free(pmom[i]);
		free(pSavEVDist[i]);
	}
	if(NUMMAT<14 && NUMMAT>1 && compC!=0)
	{
		for(i=0;i<NUMMAT;i++)
			free(gammamat[i]);
		free(pgammamult);
	}
	free(numev);
	free(SmallEV);
	free(LargEV);
	free(BinsEV);
	free(pComm);
	free(piComm);
	free(SavEVCommDist);
	free(pphi);

	return 0;
}
