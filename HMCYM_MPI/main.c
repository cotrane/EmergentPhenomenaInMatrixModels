#include <mpi.h>
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
 * parallelized version using openMP; only useful for number of matrices D >= 10; otherwise identical to HMCYM code
 * the individual matrix multiplication is not paralellized as this does not improve performance for matrices < 1000
 * different nodes always compute whole multiplication of 2 matrices and individual results are then added up in matrix
 * in master node with rank = 0
 */

struct tm tim;

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

// define MPI_operation to add two complex matrices
void CMPI_addmat (doublecomplex *a, doublecomplex *b, int *len, MPI_Datatype *type);

int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, MATRIX_SIZE=1, THERM=1, STEPS=1, STEP=1, NUMMAT=2;
	double EPS=1, LOWERACC=70, UPPERACC=90, BINWIDTH=0.01, BINWIDTHCOMM=0.01, BINWIDTHDIRAC=0.05;

	int i,k,t=0;

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

//	doublecomplex  **pmatXY, **pmatX2, **pmatY2;
//	double *trX2Y2, *trXYXY, *trComm2, *Corrratio;

	doublecomplex *pgammamult, *pgammamulttmp;
	int gammamultsize=0, numgammamult=0;
	double Smallgamma=0, Largegamma=0;
	int *SavGammaDist, BinsGamma=0;
	int gammasize=0, compC=0;
	doublecomplex **gammamat;

	doublecomplex *pComm, *piComm;
	double SmallEVComm=0, LargEVComm=0;
	int numevComm=0, compComm=1;
	int *SavEVCommDist, BinsEVComm=0;

	doublecomplex *pDirac;
	int numDirac=0, diracsize=0, compD=0;
	double SmallDirac=0, LargeDirac=0;
	int *SavDiracDist, BinsDirac=0;

	doublecomplex **tmpmat, *pmatX2, *pmatXY, *pmatY2;
	double mass=0, trX2=0, trXi2=0, trX2Y2=0, trXYXY=0, trXY=0, trX2Y=0;
	double meantrX2=0, meantrXY=0, meantrX2Y2=0, meantrXYXY=0, meantrX2Y=0, meantrXi2=0;

	int compGdd=0;
	doublecomplex *matGdd;
	double largeEVGdd=-100000, smallEVGdd=10000;
	int *saveEVGdd, numEVGdd=0, binsGdd=0;

	int start=0;

	FILE *out = NULL;
	FILE *out2 = NULL;
	FILE *out3 = NULL;
	FILE *out4 = NULL;
//	FILE *out5;
//	FILE *out6;
	char filename[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
//	char filename5[256];
//	char filename6[256];

	time_t m_startTime;

	// initialize number of processors used for simulation; value comes from command line argument
	MPI_Init(&argc, &argv);
	// initizializes MPI_communicator for all processors of number = numprocs
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	// gives each processor and unique identifier = rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// useless function; returns processor names in processor_name??
	MPI_Get_processor_name(processor_name, &namelen);

	if(argc>1)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &NUMMAT,
				&mass, &start, &compC, &compD, &compEV, &compComm);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &EPS, &LOWERACC, &UPPERACC, &THERM, &STEPS, &STEP, pathout, &EVstore, &NUMMAT, &mass,
				&start, &compC, &compD, &compEV, &compComm);
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

		sprintf(filename, "%sYM%dMMHMCLogfile-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, MATRIX_SIZE);
		if(mass==0)
		{
			BINWIDTH = 0.01;BINWIDTHCOMM=0.01;
			sprintf(filename2, "%sYM%dMMHMCEnergy-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
			sprintf(filename3, "%sHMCeps-%dx%d.txt", pathout, MATRIX_SIZE, LOOPNUMBER);
			sprintf(filename4, "%s%dMMtrX2-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
//			sprintf(filename5, "%sYM%dMMHMCEnergyS0+S1-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
//			sprintf(filename6, "%sYM%dMMlambdaMatrixEV-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
		}
		if(mass!=0)
		{
			BINWIDTH = 0.002;BINWIDTHCOMM=0.0002;
			sprintf(filename2, "%sYM%dMMHMCEnergy-%dx%d-m=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, mass);
			sprintf(filename3, "%sHMCeps-%dx%d-m=%.3f.txt", pathout, MATRIX_SIZE, LOOPNUMBER, mass);
			sprintf(filename4, "%s%dMMtrX2-%dx%d.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER);
//			sprintf(filename5, "%sYM%dMMHMCEnergyS0+S1-%dx%d-m=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, mass);
//			sprintf(filename6, "%sYM%dMMlambdaMatrixEV-%dx%d-m=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, mass);
		}

		out = fopen(filename, "a");
		out2 = fopen(filename2, "w");
		out3 = fopen(filename3, "w");
		out4 = fopen(filename4, "w");
//		out5 = fopen(filename5, "w");
//		out6 = fopen(filename6, "w");

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
	}

	N = MATRIX_SIZE*MATRIX_SIZE;

	for(i=0;i<parsnum;i++)
		parsname[i] = (char*) malloc(sizeof(char)*256);
	pmat = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pmatold = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	pmom = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
	numev = (int*) malloc(sizeof(int)*NUMMAT);
	SmallEV = (double*) malloc(sizeof(double)*NUMMAT);
	LargEV = (double*) malloc(sizeof(double)*NUMMAT);
//	pmatXY = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
//	pmatX2 = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
//	pmatY2 = (doublecomplex**) malloc(sizeof(doublecomplex*)*NUMMAT);
//	trComm2 = (double*) malloc(sizeof(double)*NUMMAT*NUMMAT);
//	Corrratio = (double*) malloc(sizeof(double)*NUMMAT*NUMMAT);
//	trXYXY = (double*) malloc(sizeof(double)*NUMMAT);
//	trX2Y2 = (double*) malloc(sizeof(double)*NUMMAT);
	pSavEVDist = (int**) malloc(sizeof(int*)*NUMMAT);
	BinsEV = (int*) malloc(sizeof(int)*NUMMAT);
	tmpmat = (doublecomplex**) malloc(sizeof(doublecomplex)*3);
	pComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	piComm = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	pmatXY = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	pmatX2 = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	pmatY2 = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
	for(i=0;i<NUMMAT;i++)
	{
		pmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmatold[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmom[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
//		pmatXY[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
//		pmatX2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
//		pmatY2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	}
	for(i=0;i<3;i++)
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);

//	memset(trComm2, 0, NUMMAT*NUMMAT*sizeof(double));
//	memset(Corrratio, 0, NUMMAT*NUMMAT*sizeof(double));
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

	if(NUMMAT<14)
	{
		gammasize = computeGammasize(NUMMAT);
		gammamultsize = MATRIX_SIZE*gammasize;
		diracsize = gammasize*N;
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
		gammamatrices(gammamat, NUMMAT);
		if(rank==0){
			printf("gammasize=%d\tgammamultsize=%d\n", gammasize, gammamultsize);
		}
	}

	if(rank==0 && start==1)
		gen_randphicplx(pmat, NUMMAT, MATRIX_SIZE);

	S = actionYM(pmat, MATRIX_SIZE, NUMMAT, mass, numprocs, rank);

	if(rank==0)
		printf("S=%f\n", S);

	for(k=1;k<LOOPNUMBER;k++)
	{
		if(rank==0)
		{
			ratio = (k*1.0)/(LOOPNUMBER*1.0) *100.0;												//prints percentage on screen
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

			if((k%100)==0)
			{
				if(acccheck < LOWERACC)											//dynamical fit of interval for fluctuation
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

			for(i=0;i<NUMMAT;i++)
			{
				memcpy(pmatold[i], pmat[i], sizeof(doublecomplex)*N);
			}

			gen_gaussmomcplx(pmom, NUMMAT, MATRIX_SIZE, 1.0);

			P = mom(pmom, NUMMAT, MATRIX_SIZE);
			H = S + P;
		}

		// broadcast variable EPS to each node
		if((k%100)==0)
		{
			MPI_Bcast(&EPS, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		for(i=0;i<STEPS;i++)
		{
			//Leapfrog

//			addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);
//			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS, C2, C3, LAMBDA);
//			addmat(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE);

			//Omelyan

			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, numprocs, rank);
			deltaYM(pmat, pmom, NUMMAT, MATRIX_SIZE, EPS/2.0, mass, numprocs, rank);
			addmat(pmat, pmom, EPS*(1-2*ksi), NUMMAT, MATRIX_SIZE, numprocs, rank);
			deltaYM(pmat, pmom, NUMMAT, MATRIX_SIZE, EPS/2.0, mass, numprocs, rank);
			addmat(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, numprocs, rank);
		}
		Send = actionYM(pmat, MATRIX_SIZE, NUMMAT, mass, numprocs, rank);
		if(rank==0)
		{
			Pend = mom(pmom, NUMMAT, MATRIX_SIZE);
			Hend = Send + Pend;

//			fprintf(out5, "%f\t%f\t%f\n",Hend, Send, Pend);

			deltaH = H - Hend;

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
		}

		if(k%STEP == 0 && k > THERM)
		{
			// these computations are not parallelized; not worth it
			if(rank==0)
			{
				expdH += exp(deltaH);

				trX2=0;trXi2=0;trXY=0;trX2Y2=0;trXYXY=0;trX2Y=0;
				for(i=0;i<NUMMAT;i++)																					//Radius Tr(X_aX^a)
				{
					diagMulti(&trXi2, pmat[i], pmat[i], MATRIX_SIZE, 1.0);
				}
				memset(pmatX2, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(pmatY2, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				memset(pmatXY, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
				Multi5(pmatX2, pmat[0], pmat[0], MATRIX_SIZE, 1.0);
				Multi5(pmatY2, pmat[1], pmat[1], MATRIX_SIZE, 1.0);
				Multi5(pmatXY, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
				diagMulti(&trX2, pmat[0], pmat[0], MATRIX_SIZE, 1.0);
				diagMulti(&trXY, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
				diagMulti(&trXYXY, pmatXY, pmatXY, MATRIX_SIZE, 1.0);
				diagMulti(&trX2Y2, pmatX2, pmatY2, MATRIX_SIZE, 1.0);
				diagMulti(&trX2Y, pmatX2, pmat[1], MATRIX_SIZE, 1.0);
				meantrX2+=trX2; meantrXi2+=trXi2; meantrXY+=trXY; meantrX2Y2+=trX2Y2; meantrXYXY+=trXYXY; meantrX2Y+=trX2Y;
				fprintf(out4, "%f\t%f\t%f\t%f\t%f\t%f\n",  trX2/MATRIX_SIZE, trXi2/MATRIX_SIZE, trXY/MATRIX_SIZE, trX2Y2/MATRIX_SIZE,
						trXYXY/MATRIX_SIZE, trX2Y/MATRIX_SIZE);
				fflush(out4);

				if(compEV!=0)
				{
					for(i=0;i<NUMMAT;i++)																								//Eigenvalues
						EVonthefly(&pSavEVDist[i], pmat[i], MATRIX_SIZE, BINWIDTH, &BinsEV[i], &LargEV[i], &SmallEV[i], &numev[i]);
				}
				if(compComm!=0)
				{
					memset(pComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
					memset(piComm, 0, sizeof(doublecomplex)*MATRIX_SIZE*MATRIX_SIZE);
					Comm(pComm, pmat[0], pmat[1], MATRIX_SIZE, 1.0);
					AddMat4(piComm, pComm, MATRIX_SIZE, 1.0);
					EVonthefly(&SavEVCommDist, piComm, MATRIX_SIZE, BINWIDTHCOMM, &BinsEVComm, &LargEVComm, &SmallEVComm, &numevComm);
				}
				if(compGdd!=0){
					computeGdd(matGdd, pmat, NUMMAT, MATRIX_SIZE, rank, numprocs);
					EVonthefly(&saveEVGdd, matGdd, NUMMAT, BINWIDTH, &binsGdd, &largeEVGdd, &smallEVGdd, &numEVGdd);
				}
			}
			// paralellized computations
			if(NUMMAT<14)
			{
				if(compC!=0)
				{
					memset(pgammamult, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
					memset(pgammamulttmp, 0, gammamultsize*gammamultsize*sizeof(doublecomplex));
					// split computation between different nodes
					for(i=rank;i<2;i+=numprocs)
						Multilambda(pgammamulttmp, pmat[i], gammamat[i], MATRIX_SIZE, gammasize);
					// add all results from individual nodes in pgammamulttmp in master node with rank = 0
					MPI_Reduce(pgammamulttmp, pgammamult, gammamultsize*gammamultsize, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
					if(rank==0)
						EVonthefly(&SavGammaDist, pgammamult, gammamultsize, BINWIDTH, &BinsGamma, &Largegamma, &Smallgamma, &numgammamult);
				}
				if(compD!=0)
				{
					memset(pDirac, 0, sizeof(doublecomplex)*diracsize*diracsize);
					genD(pDirac, pmat, gammamat, MATRIX_SIZE, gammasize, NUMMAT, numprocs, rank);
					if(rank==0)
						EVonthefly(&SavDiracDist, pDirac, diracsize, BINWIDTHDIRAC, &BinsDirac, &LargeDirac, &SmallDirac, &numDirac);
				}
			}
			if(rank==0)
			{
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
	}
	if(rank==0)
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
				SaveEVendofflight(SavDiracDist, BinsDirac, BINWIDTH, SmallDirac, numDirac);
				FileClose(EVstore);
			}
		}
		if(compGdd!=0){
			FileOpen(0, MATRIX_SIZE, LOOPNUMBER, EVstore, pathout, name5, parsnum, parsname, parsval);
			SaveEVendofflight(saveEVGdd, binsGdd, BINWIDTH, smallEVGdd, numEVGdd);
			FileClose(EVstore);
		}
		m_startTime = time(NULL);
		localtime_r(&m_startTime, &tim);
		printf("\n");
		printf("YM with %d matrices\n", NUMMAT);
		printf("loopnumber=%d\tsize=%d\teps=%lf\tloweracc=%lf\tuacc=%lf\ttherm=%d\tintegration steps=%d\tevery %d step saved\n",
				LOOPNUMBER, MATRIX_SIZE, EPS, LOWERACC, UPPERACC, THERM, STEPS, STEP);
		printf("pathout=%s\n", pathout);
		printf("acc=%d\tacc rate = %.2lf\n", acc, 100.0*acc/(k*1.0));
		printf("<expdH>=%f\t<H>=%f\t<S>=%f\t<trXi2/N>=%f\n", expdH/(therm*1.0), meanH/(therm*1.0), meanS/(therm*1.0), meantrXi2/(therm*1.));
		printf("ID = %f\n", (4*meanS)/(therm*1.0*(N-1.0)));
		printf("<trX2/N>=%f\t<trXi2/N>=%f\t<trXY/N>=%f\t<trX2Y2/N>=%f\t<trXYXY/N>=%f\t<trX2Y/N>=%f\n",
				meantrX2/(1.*MATRIX_SIZE*therm), meantrXi2/(1.*MATRIX_SIZE*therm), meantrXY/(1.*MATRIX_SIZE*therm), meantrX2Y2/(1.*MATRIX_SIZE*therm),
				meantrXYXY/(1.*MATRIX_SIZE*therm), meantrX2Y/(1.*MATRIX_SIZE*therm));
		printf("Finished at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

		out = fopen(filename, "a");
		fprintf(out,"<expdH> = %f\tAcceptance Rate = %f\n<S> = %f\t <trXi2>/N=%f\n",
				expdH/(therm*1.0), 100.0*acc/(k*1.0), meanS/(therm*1.0), meantrXi2/(therm*1.*MATRIX_SIZE));
		fprintf(out,"<trX2/N>=%f\t<trXi2/N>=%f\t<trXY/N>=%f\t<trX2Y2/N>=%f\t<trXYXY/N>=%f\t<trX2Y/N>=%f\n",
				meantrX2/(1.*MATRIX_SIZE*therm), meantrXi2/(1.*MATRIX_SIZE*therm), meantrXY/(1.*MATRIX_SIZE*therm),	meantrX2Y2/(1.*MATRIX_SIZE*therm),
				meantrXYXY/(1.*MATRIX_SIZE*therm), meantrX2Y/(1.*MATRIX_SIZE*therm));
		fprintf(out, "ID=%f\n", (4*meanS)/(t*1.0*(N-1.0)));
		fprintf(out,"Finished at %d:%d:%d\t%d/%d\n\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	//	fprintf(out3, "%f", trX2/(therm*1.*MATRIX_SIZE));
	//	for(i=0;i<NUMMAT;i++)
	//	{
	//		for(j=0;j<NUMMAT;j++)
	//		{
	//			fprintf(out3, "\t%f\t%f", trComm2[i*NUMMAT+j]/(therm*1.), Corrratio[i*NUMMAT+j]/(therm*1.));
	//		}
	//		fprintf(out3, "\n");
	//	}

		fclose(out);
		fclose(out2);
		fclose(out3);
		fclose(out4);
//		fclose(out5);
	//	fclose(out6);
	}

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
		free(pgammamulttmp);
		free(pDirac);
	}
	MPI_Finalize();
	return 0;
}

// definition of complex addition for openmpi function calls
void CMPI_addmat(doublecomplex *a, doublecomplex *b, int *len, MPI_Datatype *type)
{
	int i;

	for(i = 0; i < *len; i++)
	{
		b[i].r = a[i].r + b[i].r;
		b[i].i = a[i].i + b[i].i;
	}
}
