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
#include "model.h"
#include "MatrixMan.h"
#include "gammamatrices.h"
#include "RandomGens.h"
#include<mpi.h>

/*
 * code performs a Hybrid-Monte-Carlo simulation of Yang-Mills-Myers model for
 * 			-) 3 matrices: all files that include 3MM
 * 			-) 8 matrices: all files that include 8MM
 * 		and for modified Yang-Mills-Myers model for 8 matrices: also includes file S1.c, deltaS1.c
 * 	which model to a run the simulation for is specified in model.h
 * 	this code uses openmpi for parallelization
 */


// struct to measure length of simulation
struct tm tim;

// defines binwidth that is used for eigenvalue distributions
#define BINWIDTH	0.01
#define BINWIDTHDIRAC	0.05
#define BINWIDTHCOMM 0.01
// shortcut to print linenumber - for debugging
#define linePrint() printf("line %d\n", __LINE__)

#include <sys/types.h>
#include <sys/stat.h>

// define MPI_datatype for complex matrices
MPI_Datatype	CMPI_COMPLEX_TYPE;
// define MPI_operation to add two complex matrices
MPI_Op			CMPI_ADDCMAT_OP;
// define global variables needed for MPI: number of processors = numprocs, distributing work for node with var = rank,
// namelen, and processor_name needed for function get_processor name; no real use....
// mpiStatus = return value of mpi_function calls
int numprocs, rank, namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];
int mpiStatus;
int MATRIX_SIZE;

// define MPI_operation to add two complex matrices
void CMPI_addmat (doublecomplex *a, doublecomplex *b, int *len, MPI_Datatype *type);

int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, THERM=1, STEPS=1, STEP=1, conf=1;
	double ALPHATILDE=4, EPS=1/*, KAPPA=1*/;
	double LOWERACC=70, UPPERACC=90, LAMBDA=100;

	int i,k;

	int N, L;
	doublecomplex *pmat[NUMMAT], *pmatold[NUMMAT], *pmom[NUMMAT];
	double S=0, H=0, P=0, deltaH=0, Send=0, Hend=0, Pend=0, trYM=0, trCS=0, meanYM=0, meanCS=0, trYMold=0, trCSold=0;

	int acc=0;
	int acccheck=0;
	double mcrand=0;

	int ratio=0;											//for percentage on screen
	int percentage=0;

	double expdH=0, meanH=0, meanS=0;

	double *SmallEV, *LargEV;
	int *numev, **pSavEVDist, *BinsEV, compEV=0;
	double ksi=0.1931833;

	int EVstore, numrep=0, *rep, diag=0;

	double C2=1, C3;
	double S0=0,S1=0;
	double mommult=1.0, prop1=1, prop2=1;

	time_t m_startTime;

	doublecomplex *pgammamult, *pgammamulttmp;
	int gammamultsize=0, numgammamult=0;
	double Smallgamma=0, Largegamma=0;
	int *SavGammaDist, BinsGamma=0;
	int gammasize=0, compC=0;
	doublecomplex **gammamat;

	doublecomplex *pComm, *piComm;
	double SmallEVComm=0, LargEVComm=0;
	int numevComm=0;
	int *SavEVCommDist, BinsEVComm=0, compComm=0;

	doublecomplex *pComm2, *piComm2;
	double SmallEVComm2=1000000000000000, LargEVComm2=-1000000000000000;
	int numevComm2=0;
	int *SavEVCommDist2, BinsEVComm2=0, compComm2=0;

	doublecomplex *pDirac;
	int numDirac=0, diracsize=0;
	double SmallDirac=0, LargeDirac=0;
	int *SavDiracDist, BinsDirac=0, compD=0;

	doublecomplex **ppLieGen, *pB;
	int BSize, NumEVB=0;
	double SmallEVB=0, LargeEVB=0;
	int *pSavEVBDist, BinsEVB=0, compB=0;

	double trX2=0, meantrX2=0;
	int therm=0;

	char pathout[1000], name[256], name2[256], name3[256], name4[256], name5[256];
	int parsnum=2;
	char *parsname[2];
	double parsval[2];

	FILE *out = NULL;
	FILE *out2 = NULL;
//	FILE *out3 = NULL;
	FILE *out4 = NULL;
	FILE *out5 = NULL;
	char filename[256];
	char filename2[256];
//	char filename3[256];
	char filename4[256];
	char filename5[256];

	clock_t begin=0, end=0;
	int hours,min,sec;

	// initialize number of processors used for simulation; value comes from command line argument
	MPI_Init(&argc, &argv);
	// initizializes MPI_communicator for all processors of number = numprocs
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	// gives each processor and unique identifier = rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// useless function; returns processor names in processor_name??
	MPI_Get_processor_name(processor_name, &namelen);

	// loading parameters for simulation from file or command line
	if(argc > 4)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &ALPHATILDE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, &LAMBDA,
				pathout, &conf, &EVstore, &mommult, &numrep,&rep, &diag, &prop1, &prop2, &compC, &compD, &compEV, &compComm, &compComm2, &compB);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &ALPHATILDE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, &LAMBDA, pathout,
				&conf, &mommult, &EVstore, &numrep, &rep,&diag, &prop1, &prop2, &compC, &compD, &compEV, &compComm, &compComm2, &compB);
	}

	// defines the new datatype CMPI_COMPLEX_TYPE as 2 MPI_double variables
	mpiStatus = MPI_Type_contiguous(2, MPI_DOUBLE, &CMPI_COMPLEX_TYPE);
	if(mpiStatus != MPI_SUCCESS)
	{
		printf("type not initialized correctly!\n");
	}
	MPI_Type_commit(&CMPI_COMPLEX_TYPE);
	// define operator that connects with user-defined function CMPI_addmat
	MPI_Op_create(CMPI_addmat, 1, &CMPI_ADDCMAT_OP);

	if(rank==0)
	{
		// initializing seed for PRNG
		mt_goodseed();

		sprintf(filename, "%s%dMMHMCLogfile-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, MATRIX_SIZE);
		if(LAMBDA!=0)
			sprintf(filename2, "%s%dMMHMCEnergy-%dx%d-a=%.4lf-l=%.2f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
		else
			sprintf(filename2, "%s%dMMHMCEnergy-%dx%d-a=%.4lf.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE);
//		sprintf(filename3, "%sHMCeps-%dx%d-a=%.2lf-l=%.2f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
		if(LAMBDA!=0)
			sprintf(filename4, "%s%dMMtrX2-%dx%d-a=%.4lf-l=%.2f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
		else
			sprintf(filename4, "%s%dMMtrX2-%dx%d-a=%.4lf.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE);
		sprintf(filename5, "%s%dMMHMCYMCSS0S1-%dx%d-a=%.2lf-l=%.2f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA);
		out = fopen(filename, "a");
		out2 = fopen(filename2, "w");
//		out3 = fopen(filename3, "w");
		out4 = fopen(filename4, "w");
		out5 = fopen(filename5, "w");

		begin=clock();
		m_startTime = time(NULL);
		localtime_r(&m_startTime, &tim);
		printf("loopnumber=%d\tsize=%d\tatilde=%lf\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\tlambda=%f\tEVstore=%d\tnumrep=%d"
				"\t diag=%d\t prop=%f\tstart=%d\n",
				LOOPNUMBER, MATRIX_SIZE, ALPHATILDE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP, LAMBDA, EVstore, numrep, diag, prop1, conf);
		printf("reps: ");
		for(i=0;i<numrep;i++)
			printf("%d ", rep[i]);
		printf("\n");
		printf("pathout=%s\tnumprocs=%d\n", pathout, numprocs);
		printf("Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
		printf("%dMM code in use\n", NUMMAT);
		fprintf(out,"Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
		fprintf(out,"%dMM HMC Simulation with %d Loops, %d integration steps, saving every %d step\n "
				"%dx%d matrices, Alphatilde = %lf, Thermalization = %d loops, starting with Epsilon = %lf, Range of Acceptance Rate is %lf to %lf\n"
				"%d integration loops, lambda = %lf\n",
				NUMMAT, LOOPNUMBER, STEPS, STEP, MATRIX_SIZE, MATRIX_SIZE, ALPHATILDE, THERM, EPS, UPPERACC, LOWERACC, STEPS, LAMBDA);
		if(conf==0)
			fprintf(out,"Starting from SU(3) configuration\n");
		else if(conf==1)
			fprintf(out,"Starting from hot configuration\n");
		else if(conf==2 && diag==0)
			fprintf(out,"Starting from SU(2) configuration\n");
		else if(conf==2 && diag==1)
			fprintf(out,"Starting from SU(2)xU(1) configuration\n");
		else
			fprintf(out,"Starting from different SU(2) embedding\n");
		fprintf(out,"using Alphatilde = alpha * size^0.25! \n");
		fclose(out);
	}
	N = MATRIX_SIZE*MATRIX_SIZE;
	BSize = LIEGROUP*MATRIX_SIZE;
	gammasize = computeGammasize(NUMMAT);

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

	for(i=0;i<parsnum;i++)
		parsname[i] = (char*) malloc(sizeof(char)*256);

	pSavEVDist = (int**) malloc(sizeof(int*)*NUMMAT);
	BinsEV = (int*) malloc(sizeof(int)*NUMMAT);
	numev = (int*) malloc(sizeof(int)*NUMMAT);
	SmallEV = (double*) malloc(sizeof(double)*NUMMAT);
	LargEV = (double*) malloc(sizeof(double)*NUMMAT);
	pComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	piComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	pComm2 = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	piComm2 = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	SavEVCommDist2 = (int*) malloc(sizeof(int*));
	ppLieGen = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pB = (doublecomplex*) malloc(sizeof(doublecomplex)*BSize*BSize);
	for(i=0;i<NUMMAT;i++)
	{
		pmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmatold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmom[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		ppLieGen[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*LIEGROUP*LIEGROUP);
	}
	memset(numev, 0, sizeof(int)*NUMMAT);
	memset(SmallEV, 0, sizeof(double)*NUMMAT);
	memset(LargEV, 0, sizeof(double)*NUMMAT);
	memset(BinsEV, 0, sizeof(int)*NUMMAT);
	for(i=0;i<NUMMAT;i++)
	{
		memset(pmat[i], 0, sizeof(doublecomplex)*N);
		memset(pmatold[i], 0, sizeof(doublecomplex)*N);
		memset(pmom[i], 0, sizeof(doublecomplex)*N);
		memset(ppLieGen[i], 0, sizeof(doublecomplex)*LIEGROUP*LIEGROUP);
	}

	if(NUMMAT<14)
	{
		gammamultsize = MATRIX_SIZE*gammasize;
		diracsize = MATRIX_SIZE*MATRIX_SIZE*gammasize;
		pgammamult = (doublecomplex*) malloc(gammamultsize*gammamultsize*sizeof(doublecomplex));
		pgammamulttmp = (doublecomplex*) malloc(gammamultsize*gammamultsize*sizeof(doublecomplex));
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
		// obtain gammamatrices for correct dimensionality (d=3 or d=8)
		gammamatrices(gammamat, NUMMAT);
		printf("gammasize=%d\tgammamultsize=%d\n", gammasize, gammamultsize);
	}
	// obtain lie generators for correct dimensionality (d=3 or d=8)
	RepsSUn(ppLieGen, LIEGROUP, LIEGROUP, NUMMAT, 2);
	if(rank==0)
		printmat("lie generators:", ppLieGen[0], LIEGROUP);

	if(rank==0)
	{
		parsname[0] = "a";
		parsname[1] = "l";
		parsval[0] = ALPHATILDE;
		parsval[1] = LAMBDA;

		// define starting configuration
		startconf(pmat, ALPHATILDE, NUMMAT, MATRIX_SIZE, LIEGROUP, numrep, rep, conf, &L, diag, prop1, prop2);

		// compute Casimir operators; used for modified 8MM
		C2=(LIEGROUP-1.0)/(2.0*LIEGROUP) * L*(L+LIEGROUP);
		C3=/*C2* */((LIEGROUP*1.0-2.0)*(2.0*L+LIEGROUP*1.0))/(2.0*LIEGROUP);
	}
	// Broadcast Casimir operators to all nodes
	MPI_Bcast(&C2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&C3, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// compute value of action for current model
	S = action_function(pmat, &trYM, &trCS, ALPHATILDE, MATRIX_SIZE, NUMMAT, C2, C3, LAMBDA, &S0, &S1, numprocs, rank);										//need to change for model change!!!
	trYMold=trYM; trCSold=trCS;
	if(rank==0)
	{
		fprintf(out2,"%f\n", S);
		printf("S=%f\tS0=%f\tS1=%f\n", S,S0,S1);
	}

	// starting MC-loops
	for(k=1;k<LOOPNUMBER;k++)
	{
		if(rank==0)
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
			gen_gaussmomcplx(pmom, NUMMAT, MATRIX_SIZE, mommult);

			// computing initial momentum and hamiltonian
			P = mom(pmom, NUMMAT, MATRIX_SIZE, LAMBDA);
			H = S + P;
			trYM=0;trCS=0;
		}

		// broadcast variable EPS to each node
		if((k%100)==0)
		{
			MPI_Bcast(&EPS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		// following hamiltonian equations of motion for number of steps = STEPS
		for(i=0;i<STEPS;i++)
		{
			//Leapfrog
			//	addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);
			//	deltaX3MM(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS);
			//			addmom(pmom, pmat, EPS, NUMMAT, MATRIX_SIZE, KAPPA);
			//	addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);

			//Omelyan
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, numprocs, rank);
			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS/2.0, C2, C3, LAMBDA, numprocs, rank);
			addmat(pmat, pmom, EPS*(1-2*ksi), NUMMAT, MATRIX_SIZE, numprocs, rank);
			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS/2.0, C2, C3, LAMBDA, numprocs, rank);
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, numprocs, rank);
		}

		// computation of final value of action
		Send = action_function(pmat, &trYM, &trCS, ALPHATILDE, MATRIX_SIZE, NUMMAT, C2, C3, LAMBDA, &S0, &S1, numprocs, rank);

		if(rank==0)
		{
			// computation of final value of momentum and hamiltonian
			Pend = mom(pmom, NUMMAT, MATRIX_SIZE, LAMBDA);
			Hend = Send + Pend;

			// comparing initial and final hamiltonian
			deltaH = H - Hend;

			//  ---------------- Metropolis step ---------------------------------
			if(deltaH>0)
			{
				H = Hend;
				S = Send;
				P = Pend;
				trYMold=trYM; trCSold=trCS;
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
					trYMold=trYM; trCSold=trCS;
					acc++;
					acccheck++;
				}
				else
				{
					trYM=trYMold; trCS=trCSold;
					for(i=0;i<NUMMAT;i++)
					{
						memcpy(pmat[i], pmatold[i], sizeof(doublecomplex)*N);
					}
				}
			}
			// save different terms of action in file
			fprintf(out2,"%f\n", S);
			fprintf(out5, "%f\t%f\t%f\t%f\n",trYM, trCS, S0, S1);
		}

		// computation of observables every STEP steps if loopnumber is larger thermalization threshold
		if((k%STEP == 0) && (k > THERM))
		{
			if(rank==0)
			{
				// compute code-checker exp(\delta H); should be close to 1
				expdH += exp(deltaH);

				// compute EV dist for matrices and save in pSavEVDist or store EV's in array SavEVStore
				if(compEV!=0)
				{
					for(i=0;i<NUMMAT;i++)
						EVonthefly(&pSavEVDist[i], pmat[i], MATRIX_SIZE, BINWIDTH, &BinsEV[i], &LargEV[i], &SmallEV[i], &numev[i]);
				}

				// compute \sum_{\mu} Tr(X_{\mu}^2)
				trX2=0;
				for(i=0;i<NUMMAT;i++)
				{
					diagMulti(&trX2, pmat[i], pmat[i], MATRIX_SIZE, 1.0);
				}
				meantrX2+=trX2;
				fprintf(out4, "%f\n", trX2/MATRIX_SIZE);

				// compute EV dist for commutator i[X_1,X_2] or store EV's
				if(compComm != 0){
					memset(pComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
					memset(piComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
					Comm(pComm, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
					AddMat4(piComm, pComm, MATRIX_SIZE, 1.0);
					EVonthefly(&SavEVCommDist, piComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm, &LargEVComm, &SmallEVComm, &numevComm);
				}
				// compute EV dist for commutator i[X_1,X_8] or store EV's
				if(compComm2 != 0){
					memset(pComm2, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
					memset(piComm2, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
					Comm(pComm2, pmat[4], pmat[7], MATRIX_SIZE, 1.0);
					AddMat4(piComm2, pComm2, MATRIX_SIZE, 1.0);
					EVonthefly(&SavEVCommDist2, piComm2, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm2, &LargEVComm2, &SmallEVComm2, &numevComm2);
				}
				// compute EV distribution for matrix B = \sigma_{\mu} \otimes X_{\mu} where \mu=1,...,d
				if(compB!=0)
				{
					memset(pB, 0, sizeof(doublecomplex)*BSize*BSize);
					for(i=0;i<NUMMAT;i++)
					{
						Multilambda(pB, pmat[i], ppLieGen[i], MATRIX_SIZE, LIEGROUP);
					}
					EVonthefly(&pSavEVBDist, pB, BSize, BINWIDTH, &BinsEVB, &LargeEVB, &SmallEVB, &NumEVB);
				}
			}
			if(NUMMAT<14)
			{
				// compute EV dist of matrix C = \gamma_{\mu} \otimes X_{\mu} where \mu = 1,...,d
				if(compC != 0)
				{
					memset(pgammamult, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
					memset(pgammamulttmp, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
					for(i=rank;i<NUMMAT;i+=numprocs)
						Multilambda(pgammamulttmp, pmat[i], gammamat[i], MATRIX_SIZE, gammasize);
					MPI_Reduce(pgammamulttmp, pgammamult, gammamultsize*gammamultsize, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
					if(rank==0)
						EVonthefly(&SavGammaDist, pgammamult, gammamultsize, BINWIDTH, &BinsGamma, &Largegamma, &Smallgamma, &numgammamult);
				}
				// compute EV dist of Dirac operator D = \gamma_{\mu} \otimes [X_{\mu}, . ] where \mu = 1,..,d
				if(compD != 0)
				{
					memset(pDirac, 0, sizeof(doublecomplex)*diracsize*diracsize);
					genD(pDirac, pmat, gammamat, MATRIX_SIZE, gammasize, NUMMAT, numprocs, rank);
					if(rank==0){
						EVonthefly(&SavDiracDist, pDirac, diracsize, BINWIDTHDIRAC, &BinsDirac, &LargeDirac, &SmallDirac, &numDirac);
					}
				}
			}
			therm++;
			if(rank==0)
			{
				/* Save EV distributions every 5000th step */
				if((therm%5000)==0)
				{
					if(compEV != 0)
					{
						for(i=0;i<NUMMAT;i++)
						{
							FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
							SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
							FileClose(EVstore);
						}
					}
					if(compComm != 0)
					{
						FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
						SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
						FileClose(EVstore);
					}
					if(compComm2 != 0){
						FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
						SaveEVendofflight(SavEVCommDist2, BinsEVComm2, BINWIDTHCOMM, SmallEVComm2, numevComm2);
						FileClose(EVstore);
					}
					if(compB!=0)
					{
						FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
						SaveEVendofflight(pSavEVBDist, BinsEVB, BINWIDTH, SmallEVB, NumEVB);
						FileClose(EVstore);
					}
					if(NUMMAT<14)
					{
						if(compC != 0)
						{
							FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name2, parsnum, parsname, parsval);
							SaveEVendofflight(SavGammaDist, BinsGamma, BINWIDTH, Smallgamma, numgammamult);
							FileClose(EVstore);
						}
						if(compD != 0)
						{
							FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name4, parsnum, parsname, parsval);
							SaveEVendofflight(SavDiracDist, BinsDirac, BINWIDTHDIRAC, SmallDirac, numDirac);
							FileClose(EVstore);
						}
					}
				}
				meanH += H; meanS += S;
				meanYM += trYM; meanCS += trCS;
			}
		}
	}
	if(rank==0)
	{
		// save final EV distribution results
		printf("\n");
		if(compEV != 0)
		{
			for(i=0;i<NUMMAT;i++)
			{
				FileOpen(i, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name, parsnum, parsname, parsval);
				SaveEVendofflight(pSavEVDist[i], BinsEV[i], BINWIDTH, SmallEV[i], numev[i]);
				FileClose(EVstore);
			}
		}
		if(compComm != 0)
		{
			FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
			SaveEVendofflight(SavEVCommDist, BinsEVComm, BINWIDTHCOMM, SmallEVComm, numevComm);
			FileClose(EVstore);
		}
		if(compComm2 != 0){
			FileOpen(1, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name3, parsnum, parsname, parsval);
			SaveEVendofflight(SavEVCommDist2, BinsEVComm2, BINWIDTHCOMM, SmallEVComm2, numevComm2);
			FileClose(EVstore);
		}
		if(compB!=0)
		{
			FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
			SaveEVendofflight(pSavEVBDist, BinsEVB, BINWIDTH, SmallEVB, NumEVB);
			FileClose(EVstore);
		}
		if(NUMMAT<14)
		{
			if(compC != 0)
			{
				FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name2, parsnum, parsname, parsval);
				SaveEVendofflight(SavGammaDist, BinsGamma, BINWIDTH, Smallgamma, numgammamult);
				FileClose(EVstore);
			}
			if(compD != 0)
			{
				FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name4, parsnum, parsname, parsval);
				SaveEVendofflight(SavDiracDist, BinsDirac, BINWIDTH, SmallDirac, numDirac);
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
		printf("ID = %f\t <trX2>/N=%f\n", (4*meanYM + 3*meanCS)/(therm*1.0*(N-1.0)), meantrX2/(therm*1.*MATRIX_SIZE));
		printf("Finished at %d:%d:%d\t%d/%d\nTime elapsed: %d h:%d min:%d s \n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1,hours, min,sec);

		out = fopen(filename, "a");
		fprintf(out,"<expdH> = %f\tAcceptance Rate = %f\n<S> = %f\t <trX2>/N=%f\n",
				expdH/(therm*1.0), 100.0*acc/(k*1.0), meanS/(therm*1.0), meantrX2/(therm*1.*MATRIX_SIZE));
		if(LIEGROUP==2)
			fprintf(out, "SU(2) ground state = %f\n", -N*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12);
		if((LIEGROUP-2)!=0)
		{
			fprintf(out, "SU(3) ground state = %f\tSU(2) ground state = %f\n",
					-MATRIX_SIZE*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12,
					-MATRIX_SIZE*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*(N-1)/24);
		}
		fprintf(out, "ID=%f\n", (4*meanYM + 3*meanCS)/(therm*1.0*(N-1.0)));
		fprintf(out,"Finished at %d:%d:%d\t%d/%d\nTime elapsed: %d h:%d min:%d s \n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1,hours, min,sec);

		fclose(out);
		fclose(out2);
		fclose(out4);
		fclose(out5);
	}

	for(i=0;i<NUMMAT;i++)
	{
//		free(pSavEVDist[i]);
		free(pmat[i]);
		free(pmatold[i]);
		free(pmom[i]);
	}
//	free(pSavEVDist);
	free(BinsEV);
	free(numev);
	free(SmallEV);
	free(LargEV);
	free(pComm);
	free(piComm);
	free(pComm2);
	free(piComm2);
	for(i=0;i<NUMMAT;i++)
		free(ppLieGen[i]);
	if(NUMMAT<14)
	{
		free(pgammamult);
		free(pgammamulttmp);
		for(i=0;i<NUMMAT;i++)
		{
			free(gammamat[i]);
		}
		if(compD!=0)
		{
			free(pDirac);
		}
	}
	free(ppLieGen);
	free(pB);

	MPI_Finalize();
	return 0;
}

// definition of complex addition for openmpi function calls
void CMPI_addmat(doublecomplex *a, doublecomplex *b, int *len, MPI_Datatype *type)
{
	int i;

	for(i = 0; i < *len; i++)
	{
		//		b[i] = a[i] + b[i];
		b[i].r = a[i].r + b[i].r;
		b[i].i = a[i].i + b[i].i;
	}
}
