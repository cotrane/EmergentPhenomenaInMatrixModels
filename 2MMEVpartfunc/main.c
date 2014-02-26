#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
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
#include <sys/types.h>
#include <sys/stat.h>

struct tm tim;

FILE *outext;
FILE *outTraces;
char filenameext[256];

#define BINWIDTHDIRAC	0.05

/*
 * code for HMC simulation of the 2-matrix model:
 * 		1) S[mass,pmat] = MATRIX_SIZE * Tr( pmat_{\mu}^2 + \frac{1}{2*mass*mass} [pmat_{\mu},pmat_{\nu}]^2 )		or
 * 		2) S[g, pmat] = MATRIX_SIZE * Tr( pmat_{\mu}^2 + g^2 [pmat_{\mu},pmat_{\nu}]^2 )
 * where \mu,\nu = 1,2.
 *
 * matrices pmat_{\mu} are diagonalized and one is integrated out to increase speed of simulation. one can re-generate the
 * second matrix as gaussian, which is done in the end (see my thesis)
 *
 * set either parameter 'mass' or 'g' unequal zero otherwise code will not work!
 *
 * code works for matrices with or without trace; current version works for matrices with trace; for traceless matrices
 * one needs to change function calls in HMC routine
 */

int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, MATRIX_SIZE=1, THERM=1, STEPS=1, STEP=1, NUMMAT=2;
	double EPS=1, LOWERACC=70, UPPERACC=90, BINWIDTH=0.002, BINWIDTHCOMM=0.0002;

	int i,j,k,l,m,t=0;

	double *pmat, *pmatold, *pmom;
	double S=0, H=0, P=0, deltaH=0, Send=0, Hend=0, Pend=0;

	int acc=0;
	int acccheck=0;
	double mcrand=0;

	int ratio=0;											//for percentage on screen
	int percentage=0;

	double expdH=0, meanH=0, meanS=0;
	double ksi=0.1931833;

	char pathout[1000], name[256], name2[256], name3[256], name4[256], name5[256], name6[256];
	char name7[256], name8[256], name9[256], name10[256], name14[256];
	int parsnum=0;
	char *parsname[1];
	double parsval[1];
	int EVstore=0;

	doublecomplex **pmatXYZ;

	int therm=0;

	double mass=0, g=0, trX2=0, trX2mean=0;
	int compTraces=0;

	doublecomplex *pgammamult;
	int gammamultsize, numgammamult=0;
	double Smallgamma=0, Largegamma=0;
	int *SavGammaDist, BinsGamma=0;
	int gammasize=2, compC=1;
	doublecomplex **gammamat;

	double *SmallEVXYZ, *LargEVXYZ, *pevX;
	int *numevXYZ, **pSavEVDist, *BinsEV;

//	doublecomplex *pmattmp, *pmattmp2;										//for testing purposes
//	double LargEVtmp=0, SmallEVtmp=0, LargEVtmp2=0, SmallEVtmp2=0;
//	double *pSavEVtmp, *pSavEVtmp2;
//	int numevtmp=0, numevtmp2=0;

	doublecomplex *pComm, *piComm;
	double SmallEVComm=0, LargEVComm=0;
	int numevComm=0;
	int *SavEVCommDist, BinsEVComm=0, compComm=0;

	doublecomplex *pAComm;
	double SmallEVAComm=0, LargEVAComm=0;
	int numevAComm=0;
	int *SavEVACommDist, BinsEVAComm=0, compAComm=0;

	doublecomplex *pDirac, *pDiracMinus, **pAdmat;
	int numDirac=0, diracsize=0;
	double SmallDirac=0, LargeDirac=0;
	int *SavDiracDist, BinsDirac=0;
	int numbasis, compD=0;

	doublecomplex *pphi;
	int *pSavDistnHr, *pSavDistnHi, BinsEVnH[2];
	int numEVnH=0, compPhi=0;
	doublecomplex LargEVnH, SmallEVnH;
	int *pSavDistMod, BinsEVMod=0, numEVMod=0;
	double LargeEVMod=0, SmallEVMod=0;

	double *SavOffDiagY;
	double largeOffDiagY, smallOffDiagY=0;
	int numOffDiagY=0;

	double *SavOffDiagY1Q;
	double largeOffDiagY1Q, smallOffDiagY1Q=0;
	int numOffDiagY1Q=0, compOffDiag=0;

	double *SavOffDiagY3Q;
	double largeOffDiagY3Q, smallOffDiagY3Q=0;
	int numOffDiagY3Q=0;

	double mult=1.;

	FILE *out;
	FILE *out2;
	FILE *out3;
//	FILE *out4;
	char filename[256];
	char filename2[256];
	char filename3[256];
//	char filename4[256];
	char filenameTr[256];

	/*
	 * read in parameters from input file or command-line (as used in python file)
	 */
	if(argc>1)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &mass, &g, &mult,
				&compComm, &compC, &compD, &compTraces, &compPhi, &compAComm, &compOffDiag, &NUMMAT);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &mass, &g, &mult,
				&compComm, &compC, &compD, &compTraces, &compPhi, &compAComm, &compOffDiag, &NUMMAT);
	}

	sprintf(name, "%d", NUMMAT);
	strncat(name, "MMEV", 4);
	sprintf(name2, "%d", NUMMAT);
	strncat(name2, "MMCEV", 5);
	sprintf(name3, "%d", NUMMAT);
	strncat(name3, "MMCommEVXY", 10);
	sprintf(name4, "%d", NUMMAT);
	strncat(name4, "MMDiracEV", 9);
	sprintf(name5, "%d", NUMMAT);
	strncat(name5, "MMphiEV", 7);
	sprintf(name6, "%d", NUMMAT);
	strncat(name6, "MMmodEV", 7);
	sprintf(name7, "%d", NUMMAT);
	strncat(name7, "MMoffDiagY", 10);
	sprintf(name8, "%d", NUMMAT);
	strncat(name8, "MMoffDiagY1Q", 12);
	sprintf(name9, "%d", NUMMAT);
	strncat(name9, "MMoffDiagY3Q", 12);
	sprintf(name10, "%d", NUMMAT);
	strncat(name10, "MMoffDiagComm", 13);
	sprintf(name14, "%d", NUMMAT);
	strncat(name14, "MMACommEV", 9);

	if(g!=0 && mass!=0)
	{
		printf("mass!=0 and g!=0 - i don't know what to do!! decide and set one to zero!\n");
		return 0;
	}

	sprintf(filename, "%s2MMLogfile-%dx%d.txt", pathout, MATRIX_SIZE, MATRIX_SIZE);
	if(g!=0)
		sprintf(filename2, "%s2MMEnergy-%dx%d-g=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, g);
	if(mass!=0)
		sprintf(filename2, "%s2MMEnergy-%dx%d-m=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, mass);
	if(g!=0)
		sprintf(filename3, "%s2MMeps-%dx%d-g=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, g);
	if(mass!=0)
		sprintf(filename3, "%s2MMeps-%dx%d-m=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, mass);
//	sprintf(filename4, "%sYM2MMtraces-%dx%d-m=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, mass);
	out = fopen(filename, "a");
	out2 = fopen(filename2, "w");
	out3 = fopen(filename3, "w");
//	out4 = fopen(filename4, "w");
	if(compTraces!=0)
	{
		sprintf(filenameTr, "%s2MMtraces2468-%dx%d-m=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, mass);
		outTraces = fopen(filenameTr, "w");
	}

	// start timing
	time_t m_startTime;
	m_startTime = time(NULL);

	// initializing seed for PRNG
	mt_goodseed();

	localtime_r(&m_startTime, &tim);
	printf("Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	printf("effective 2MM simulation with %d matrices\n", NUMMAT);
	printf("g=%f\tmass=%f\tloopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\tmult=%f\n",
			g,mass, LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP, mult);
	printf("pathout=%s\n", pathout);

	fprintf(out,"Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	fprintf(out,"YM HMC Simulation with %d matrices for %d loops, %d integration steps, saving every %d step\n "
			"%dx%d matrices, Thermalization = %d loops, starting with Epsilon = %lf, Range of Acceptance Rate is %lf to %lf\n"
			"%d integration loops\n mass=%f g=%f\n",
			NUMMAT, LOOPNUMBER, STEPS, STEP, MATRIX_SIZE, MATRIX_SIZE, THERM, EPS, UPPERACC, LOWERACC, STEPS, mass, g);
	fclose(out);

	gammamultsize = MATRIX_SIZE*gammasize;
	numbasis = MATRIX_SIZE*MATRIX_SIZE-1;
	diracsize = MATRIX_SIZE*MATRIX_SIZE*gammasize;
	LargEVnH.r=0;LargEVnH.i=0;SmallEVnH.r=0;SmallEVnH.i=0;
	pmat = (double*) malloc(sizeof(double)*MATRIX_SIZE);
	memset(pmat, 0, sizeof(double)*MATRIX_SIZE);
	pmatold = (double*) malloc(sizeof(double)*MATRIX_SIZE);
	memset(pmatold, 0, sizeof(double)*MATRIX_SIZE);
	pmom = (double*) malloc(sizeof(double)*MATRIX_SIZE);
	memset(pmom, 0, sizeof(double)*MATRIX_SIZE);
	pgammamult = (doublecomplex*) malloc(gammamultsize*gammamultsize*sizeof(doublecomplex));
	gammamat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pmatXYZ = (doublecomplex**) malloc(sizeof(doublecomplex)*NUMMAT);
	pevX = (double*) malloc(sizeof(double)*MATRIX_SIZE);
	numevXYZ = (int*) malloc(sizeof(int)*NUMMAT);
	memset(numevXYZ, 0, sizeof(int)*NUMMAT);
	SmallEVXYZ = (double*) malloc(sizeof(double)*NUMMAT);
	memset(SmallEVXYZ, 0, sizeof(double)*NUMMAT);
	LargEVXYZ = (double*) malloc(sizeof(double)*NUMMAT);
	memset(LargEVXYZ, 0, sizeof(double)*NUMMAT);
	pSavEVDist = (int**) malloc(sizeof(int*)*NUMMAT);
	pComm = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	pAComm = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	piComm = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	BinsEV = (int*) malloc(sizeof(int)*NUMMAT);
	memset(BinsEV, 0, sizeof(int)*NUMMAT);
	pDirac = (doublecomplex*) malloc(diracsize*diracsize*sizeof(doublecomplex));
	pDiracMinus = (doublecomplex*) malloc(diracsize*diracsize*sizeof(doublecomplex));
	pAdmat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pphi = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	SavOffDiagY = (double*) malloc(sizeof(double)*MATRIX_SIZE*5000);
	memset(SavOffDiagY, 0, sizeof(double)*MATRIX_SIZE*5000);
	SavOffDiagY1Q = (double*) malloc(sizeof(double)*MATRIX_SIZE*5000);
	memset(SavOffDiagY1Q, 0, sizeof(double)*MATRIX_SIZE*5000);
	SavOffDiagY3Q = (double*) malloc(sizeof(double)*MATRIX_SIZE*5000);
	memset(SavOffDiagY3Q, 0, sizeof(double)*MATRIX_SIZE*5000);
	for(i=0;i<NUMMAT;i++)
	{
		pmatXYZ[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
		memset(pmatXYZ[i], 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
		gammamat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*gammasize*gammasize);
		memset(gammamat[i], 0, sizeof(doublecomplex)*gammasize*gammasize);
		pAdmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*numbasis*numbasis);
	}
	memset(BinsEVnH, 0, sizeof(int)*2);
	memset(pmatXYZ[0], 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);

	// generate \gamma-matrices for computation of observables in the end
	gammamatrices(gammamat, NUMMAT);

//	for(i=0;i<NUMMAT;i++)	printmat("gammamat:", gammamat[i], computeGammasize(NUMMAT));

	// generate random diagonal matrix pmat as initial configuration
	gen_randphicplxtrace(pmat, MATRIX_SIZE, mult);

	// compute value of either S[mass, pmat] or S[g, pmat]
	if(mass!=0)
		S = actionYM(pmat, MATRIX_SIZE, mass);
	if(g!=0)
		S = actionYM_N(pmat, MATRIX_SIZE, g);

	fprintf(out2,"%f\n", S);
	printf("S=%f\n", S);

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

		// dynamical fit of interval for fluctuation; increases/reduces step size EPS in hamiltonian eom; done every 100th mc-step
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
			fprintf(out3, "%f\n", EPS);
		}

		// copy current matrix pmat to pmatold
		memcpy(pmatold, pmat, sizeof(double)*MATRIX_SIZE);

		// generate gaussian diagonal matrix pmom as momentum in HMC routine
		gen_gaussmomcplxtrace(pmom, MATRIX_SIZE, 1.0);

		// compute value of Tr(0.5 * pmom^2)
		P = momtrace(pmom, MATRIX_SIZE);
		// compute initial hamiltonian
		H = S + P;

		/*
		 * follow hamiltionian equations of motion to obtain final configuration of pmat which is proposed in
		 * metropolis step
		 * ksi is constant definied in Omelyan integrator
		 * STEPS is number of integration steps
		 *
		*/
		for(i=0;i<STEPS;i++)
		{
			addmattrace(pmat, pmom, ksi*EPS, MATRIX_SIZE);
			if(mass!=0)
				deltaYM(pmat, pmom, MATRIX_SIZE, EPS/2.0, mass);
			if(g!=0)
				deltaYM_N(pmat, pmom, MATRIX_SIZE, EPS/2.0, g);
			addmattrace(pmat, pmom, EPS*(1-2*ksi), MATRIX_SIZE);
			if(mass!=0)
				deltaYM(pmat, pmom, MATRIX_SIZE, EPS/2.0, mass);
			if(g!=0)
				deltaYM_N(pmat, pmom, MATRIX_SIZE, EPS/2.0, g);
			addmattrace(pmat, pmom, ksi*EPS, MATRIX_SIZE);
		}
		// compute final value of S[mass,pmat] or S[g, pmat]
		if(mass!=0)
			Send = actionYM(pmat, MATRIX_SIZE, mass);
		if(g!=0)
			Send = actionYM_N(pmat, MATRIX_SIZE, g);

		// compute final version of Tr(0.5 * pmom^2)
		Pend = momtrace(pmom, MATRIX_SIZE);
		// compute value of final hamiltonian
		Hend = Send + Pend;
		// compute difference between initial and final hamiltonian
		deltaH = H - Hend;

		// ----------------- Metropolis step: ----------------------------

		// accept if Hend < H
		if(deltaH>0)
		{
			H = Hend;
			S = Send;
			P = Pend;
			acc++;
			acccheck++;
		}
		// compare with random number in case Hend > H
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
				memcpy(pmat, pmatold, sizeof(double)*MATRIX_SIZE);
			}
		}

		fprintf(out2,"%f\n", S);
		fflush(out2);

		// compute Observables every STEP mc-steps and if current mc-step 'k' is larger a thermalization threshold THERM
		if(k%STEP == 0 && k > THERM)
		{
			// compute \sum e^{\Delta H}; should be \sim 1 when averaged over all mc-steps (averaging is done in the end)
			expdH += exp(deltaH);

			// copy pmat into pmatXYZ[0] and generate second matrix pmatXYZ[1]
			memset(pmatXYZ[0], 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
			for(i=0;i<MATRIX_SIZE;i++)
			{
				(pmatXYZ[0]+i*MATRIX_SIZE+i)->r = pmat[i];
			}
			findEigenvalues(pevX, pmatXYZ[0], MATRIX_SIZE);
			for(i=0;i<MATRIX_SIZE;i++)
			{
				(pmatXYZ[0]+i*MATRIX_SIZE+i)->r = pevX[i];
			}
			GenGaussMatrixY(pmatXYZ[1], pevX, mass, MATRIX_SIZE);

			// compute EV's for the two matrices and store their distribution
			for(i=0;i<NUMMAT;i++)
				EVonthefly(&pSavEVDist[i], pmatXYZ[i], MATRIX_SIZE, BINWIDTH, &BinsEV[i], &LargEVXYZ[i], &SmallEVXYZ[i], &numevXYZ[i]);

			// compute matrix pphi = pmatXYZ[0] + i*pmatXYZ[1] and EV distribution of real and imaginary part as well as of modulus
			if(compPhi!=0){
				memset(pphi, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				AddMat9(pphi, pmatXYZ[0], pmatXYZ[1], MATRIX_SIZE, 1.0);
				EVtoDist_nonHermitian(&pSavDistnHr, &pSavDistnHi, pphi, MATRIX_SIZE, BINWIDTH, BinsEVnH, &LargEVnH,
						&SmallEVnH, &numEVnH);
				EVtoDist_Modulus_nH(&pSavDistMod, pphi, MATRIX_SIZE, BINWIDTH, &BinsEVMod, &LargeEVMod, &SmallEVMod,
						&numEVMod);
			}
			// compute Anticommutator {pmatXYZ[0],pmatXYZ[1]} and its EV distribution
			if(compAComm!=0)
			{
				memset(pAComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
//				Multi5(pAComm, pmatXYZ[0], pmatXYZ[1], MATRIX_SIZE, 1.0);
//				Multi5(pAComm, pmatXYZ[1], pmatXYZ[0], MATRIX_SIZE, 1.0);
				AntiComm(pAComm, pmatXYZ[0], pmatXYZ[1], MATRIX_SIZE, 1.0);
				EVonthefly(&SavEVACommDist, pAComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVAComm, &LargEVAComm, &SmallEVAComm, &numevAComm);

			}
			// compute matrix C = \gamma_{\mu} \otimes pmatXYZ_{\mu}; \mu=1,2; and its EV distribution
			if(compC!=0)
			{
				memset(pgammamult, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
				for(i=0;i<NUMMAT;i++)
				{
					Multilambda(pgammamult, pmatXYZ[i], gammamat[i], MATRIX_SIZE, gammasize);
				}
				EVonthefly(&SavGammaDist, pgammamult, gammamultsize, BINWIDTH, &BinsGamma, &Largegamma, &Smallgamma, &numgammamult);
			}
			// compute i[pmatXYZ[0],pmatXYZ[1]] and its EV distribution
			if(compComm!=0)
			{
				memset(pComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(piComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				Comm(pComm, pmatXYZ[0], pmatXYZ[1], MATRIX_SIZE, 1.0);
				AddMat4(piComm, pComm, MATRIX_SIZE, 1.0);
				EVonthefly(&SavEVCommDist, pComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm, &LargEVComm, &SmallEVComm, &numevComm);
			}

			/*
			 * save modulus of off-diagonal of matrix Y to check for commuting modes in the center of the matrix
			 * no inspiring result so i didn't really use it
			*/
			if(compOffDiag!=0){
				/* main off-diagonal */
				for(i=0; i<MATRIX_SIZE;i++){
					SavOffDiagY[numOffDiagY+i] = sqrt( pow((pmatXYZ[1]+(i*MATRIX_SIZE + MATRIX_SIZE-1-i))->r,2) +
							pow((pmatXYZ[1]+(i*MATRIX_SIZE + MATRIX_SIZE-1-i))->i,2) );
					if(SavOffDiagY[numOffDiagY+i] > largeOffDiagY){
						largeOffDiagY = SavOffDiagY[numOffDiagY+i];
					}
				}
				numOffDiagY += MATRIX_SIZE;

				/* off-diagonal 1/4 above/below the main off-diagonal */
				for(i=0; i<(MATRIX_SIZE/2);i++){
					SavOffDiagY1Q[numOffDiagY1Q+i] = sqrt( pow((pmatXYZ[1]+(i*MATRIX_SIZE + (MATRIX_SIZE - 1)/2-i))->r,2) +
							pow((pmatXYZ[1]+(i*MATRIX_SIZE + (MATRIX_SIZE-1)/2 - i))->i,2) );
					if(SavOffDiagY1Q[numOffDiagY1Q+i] > largeOffDiagY1Q){
						largeOffDiagY1Q = SavOffDiagY1Q[numOffDiagY1Q+i];
					}

					SavOffDiagY3Q[numOffDiagY3Q+i] = sqrt( pow((pmatXYZ[1]+(i*MATRIX_SIZE + (MATRIX_SIZE - 1)-i) + MATRIX_SIZE * (MATRIX_SIZE - 1)/2)->r,2) +
							pow((pmatXYZ[1]+(i*MATRIX_SIZE + (MATRIX_SIZE-1) - i) + MATRIX_SIZE * (MATRIX_SIZE - 1)/2)->i,2) );
					if(SavOffDiagY3Q[numOffDiagY3Q+i] > largeOffDiagY3Q){
						largeOffDiagY3Q = SavOffDiagY3Q[numOffDiagY3Q+i];
					}
				}
				numOffDiagY1Q += MATRIX_SIZE/2;
				numOffDiagY3Q += MATRIX_SIZE/2;
			}

			/*
			 * compute dirac operator i*\gamma \otimes [pmatXYZ_{\mu}, . ] where \mu=1,2 and its EV dist
			 * commented version is correct as well but much slower
			 */
			if(compD!=0)
			{
//				memset(pDiracMinus, 0, diracsize*diracsize*sizeof(doublecomplex));
//				memset(pDirac, 0, diracsize*diracsize*sizeof(doublecomplex));
//				for(i=0;i<NUMMAT;i++)
//					memset(pAdmat[i], 0, sizeof(doublecomplex)*numbasis*numbasis);
//				ConstructAdjointAction(pAdmat, pmatXYZ, MATRIX_SIZE, NUMMAT);
//				for(i=0;i<NUMMAT;i++)
//					Multilambda2(pDirac, pAdmat[i], gammamat[i], numbasis, gammasize, 1.0);
//				AddMat4(pDiracMinus, pDirac, diracsize, 1.0);
//				EVonthefly(&SavDiracDist, pDirac, diracsize, BINWIDTH, &BinsDirac, &LargeDirac, &SmallDirac, &numDirac);

				memset(pDirac, 0, sizeof(doublecomplex)*diracsize*diracsize);
				genD(pDirac, pmatXYZ, gammamat, MATRIX_SIZE, gammasize, NUMMAT);
				EVonthefly(&SavDiracDist, pDirac, diracsize, BINWIDTHDIRAC, &BinsDirac, &LargeDirac, &SmallDirac, &numDirac);
			}

			// compute <Tr(pmatXYZ[0]^2)>; averaging is done in the end below
			trX2=0;
			for(i=0;i<MATRIX_SIZE;i++)
			{
				trX2 += pmat[i]*pmat[i];
			}

			trX2mean+=trX2;
			if(compTraces!=0)
				Traces(pmatXYZ, MATRIX_SIZE, NUMMAT, LOOPNUMBER, pathout);

			therm++;

			/*
			 * save observables every 5000 steps to avoid loss of information in case of unexpected occurrences...
			 * and to save working memory space
			 */
			if((therm%5000)==0)
			{
				if(compPhi!=0){
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
					SaveEVendofflight(pSavDistnHr, BinsEVnH[0], BINWIDTH, SmallEVnH.r, numEVnH);
					FileClose(EVstore);

					FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
					SaveEVendofflight(pSavDistnHi, BinsEVnH[1], BINWIDTH, SmallEVnH.i, numEVnH);
					FileClose(EVstore);

					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name6, parsnum, parsname, parsval);
					SaveEVendofflight(pSavDistMod, BinsEVMod, BINWIDTH, SmallEVMod, numEVMod);
					FileClose(EVstore);
				}
				for(i=0;i<NUMMAT;i++)
				{
					FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
					SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEVXYZ[i], numevXYZ[i]);
					FileClose(EVstore);
				}
				if(compComm!=0)
				{
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
					SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
					FileClose(EVstore);
				}
				if(compAComm!=0)
				{
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name14, parsnum, parsname, parsval);
					SaveEVendofflight(SavEVACommDist, BinsEVAComm, BINWIDTHCOMM, SmallEVAComm, numevAComm);
					FileClose(EVstore);
				}
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
				if(compOffDiag!=0){
					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name7, parsnum, parsname, parsval);
					DataDistribution(SavOffDiagY, largeOffDiagY, smallOffDiagY, BINWIDTH, 0, numOffDiagY);
					FileClose(EVstore);
					memset(SavOffDiagY, 0, sizeof(double)*MATRIX_SIZE*5000);
					numOffDiagY = 0;
					largeOffDiagY = -1;

					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name8, parsnum, parsname, parsval);
					DataDistribution(SavOffDiagY1Q, largeOffDiagY1Q, smallOffDiagY1Q, BINWIDTH, 0, numOffDiagY1Q);
					FileClose(EVstore);
					memset(SavOffDiagY1Q, 0, sizeof(double)*MATRIX_SIZE*5000);
					numOffDiagY1Q = 0;
					largeOffDiagY1Q = -1;

					FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name9, parsnum, parsname, parsval);
					DataDistribution(SavOffDiagY3Q, largeOffDiagY3Q, smallOffDiagY3Q, BINWIDTH, 0, numOffDiagY3Q);
					FileClose(EVstore);
					memset(SavOffDiagY3Q, 0, sizeof(double)*MATRIX_SIZE*5000);
					numOffDiagY3Q= 0;
					largeOffDiagY3Q = -1;
				}
			}
			meanH += H; meanS += S;
		}
	}
	printf("100\n");

	// save observables in the end
	if(compPhi!=0){
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
		SaveEVendofflight(pSavDistnHr, BinsEVnH[0], BINWIDTH, SmallEVnH.r, numEVnH);
		FileClose(EVstore);

		FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
		SaveEVendofflight(pSavDistnHi, BinsEVnH[1], BINWIDTH, SmallEVnH.i, numEVnH);
		FileClose(EVstore);

		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name6, parsnum, parsname, parsval);
		SaveEVendofflight(pSavDistMod, BinsEVMod, BINWIDTH, SmallEVMod, numEVMod);
		FileClose(EVstore);
	}
	for(i=0;i<NUMMAT;i++)
	{
		FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
		SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEVXYZ[i], numevXYZ[i]);
		FileClose(EVstore);
	}
	if(compComm!=0)
	{
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
		SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
		FileClose(EVstore);
	}
	if(compAComm!=0)
	{
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name14, parsnum, parsname, parsval);
		SaveEVendofflight(SavEVACommDist, BinsEVAComm, BINWIDTHCOMM, SmallEVAComm, numevAComm);
		FileClose(EVstore);
	}
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
	if(compOffDiag!=0){
		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name7, parsnum, parsname, parsval);
		DataDistribution(SavOffDiagY, largeOffDiagY, smallOffDiagY, BINWIDTH, 0, numOffDiagY);
		FileClose(EVstore);

		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name8, parsnum, parsname, parsval);
		DataDistribution(SavOffDiagY1Q, largeOffDiagY1Q, smallOffDiagY1Q, BINWIDTH, 0, numOffDiagY1Q);
		FileClose(EVstore);

		FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name9, parsnum, parsname, parsval);
		DataDistribution(SavOffDiagY3Q, largeOffDiagY3Q, smallOffDiagY3Q, BINWIDTH, 0, numOffDiagY3Q);
		FileClose(EVstore);
	}

	// print some info to screen
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("YM with 2 matrices\n");
	printf("loopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\n",
			LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP);
	printf("pathout=%s\n", pathout);
	printf("<expdH>=%f\t<H>=%f\t<S>=%f\t<trX2>/N=%f\n", expdH/(therm*1.0), meanH/(therm*1.0), meanS/(therm*1.0), trX2mean/(1.*therm*MATRIX_SIZE));
	printf("acc=%d\tacc rate = %.2lf\n", acc, 100.0*acc/(k*1.0));
	printf("Finished at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	// save some info in log-file 'out'
	out = fopen(filename, "a");
	fprintf(out,"<expdH> = %f\tAcceptance Rate = %f\t<S> = %f\t<trX2>/N = %f\n",
			expdH/(therm*1.0), 100.0*acc/(k*1.0), meanS/(therm*1.0), trX2mean/(1.*therm*MATRIX_SIZE));
	fprintf(out,"Finished at %d:%d:%d\t%d/%d\n\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	fclose(out);
	fclose(out2);
	fclose(out3);
//	fclose(out4);
//	fclose(outext);

	free(pmat);
	free(pgammamult);
	free(SavGammaDist);
	free(SmallEVXYZ);
	free(LargEVXYZ);
	free(BinsEV);
	free(numevXYZ);
	free(pComm);
	free(piComm);
	for(i=0;i<NUMMAT;i++)
	{
		free(pmatXYZ[i]);
		free(gammamat[i]);
		free(pSavEVDist[i]);
		free(pAdmat[i]);
	}
	free(gammamat);
	free(pAdmat);
	free(pDirac);
	free(pDiracMinus);
	free(pSavEVDist);
	free(SavEVCommDist);
	free(SavDiracDist);
	free(pSavDistnHr);
	free(pSavDistnHi);
	free(pSavDistMod);
	free(pphi);

	return 0;
}
