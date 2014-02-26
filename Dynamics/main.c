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
#include "RandomGens.h"

struct tm tim;
#define BINWIDTH	0.01

#include <sys/types.h>
#include <sys/stat.h>

/*
 * code follows classical hamiltonian equations of motion for various models; i.e.: there is no Metropolis step;
 * otherwise the code is equal to the MC simulation
 *
 * it's best to start the code in a local/global minimum for the 3MM and 8MM, i.e. SU(2) or SU(3) symmetric state.
 * otherwise the initial momentum is too large to see any structure; that's why it doesn't really work well with
 * lambda != 0 in the SU(2) configuration; that's no minimum!
 *
 * to change between the different models one needs to change the #define parameter in model.h:
 * for 1MM (corresponding to the harmonic oscillator): MM1; H = 0.5 * pmom^2 + 0.5 * pmat^2
 * for 3MM (Yang-Mills-Myers model): MM3
 * for 8MM (YM-Myers model plus Anticommutator term): MM8
 *
 * epsilon for \alpha=0 no problem in any model; but:
 * for the 1MM epsilon is (almost) arbitrary
 * for the 3MM epsilon ~ 0.005
 * for the 8MM epsilon ~ 0.005
 */


int main(int argc, char *argv[])
{
	int LOOPNUMBER=1, MATRIX_SIZE=1, THERM=1;
	double ALPHATILDE=4, EPS=1;
	double LOWERACC=70, UPPERACC=90, LAMBDA=100;

	int i,k;

	int N;
	doublecomplex *pmat[NUMMAT], *pmom[NUMMAT];
	double S=0, H=0, P=0, Send=0, Hend=0, Pend=0;

	int ratio=0;											//for percentage on screen
	int percentage=0;

	double ksi=0.1931833;
	int L;

	double C2=1, C3;
	double S0=0,S1=0;

	char pathout[1000];
	int conf=1;
	double mommult, mass=1.0;

	FILE *out5;
	char filename5[256];

	// initializing seed for PRNG
//	mt_goodseed();

	/*
	 * read parameters from input file
	 */
	if(argc>1)
	{
		ReadInputCL(argc, argv, &LOOPNUMBER, &MATRIX_SIZE, &ALPHATILDE, &EPS, &LOWERACC, &UPPERACC, &THERM,
				&LAMBDA, pathout, &conf, &mommult, &mass);
	}
	else
	{
		ReadInputFile(&LOOPNUMBER, &MATRIX_SIZE, &ALPHATILDE, &EPS, &LOWERACC, &UPPERACC, &THERM, &LAMBDA,
				pathout, &conf, &mommult, &mass);
	}

	sprintf(filename5, "%s%dMMHMCEnergyS0+S1-%dx%d-a=%.2lf-l=%.2f-mom=%.3f.txt", pathout, NUMMAT, MATRIX_SIZE, LOOPNUMBER, ALPHATILDE, LAMBDA, mommult);
	out5 = fopen(filename5, "w");

	time_t m_startTime;
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("Started at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);
	printf("%dMM code in use\n", NUMMAT);
	printf("loopnumber=%d\tsize=%d\tatilde=%lf\teps=%lf\tloweracc=%lf\tuacc=%lf\t therm=%d\tlambda=%f start=%d\n",
			LOOPNUMBER, MATRIX_SIZE, ALPHATILDE, EPS, LOWERACC, UPPERACC, THERM, LAMBDA, conf);
	printf("pathout=%s\n", pathout);

	N = MATRIX_SIZE*MATRIX_SIZE;

	for(i=0;i<NUMMAT;i++)
	{
		pmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		pmom[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
	}

	for(i=0;i<NUMMAT;i++)
	{
		memset(pmat[i], 0, sizeof(doublecomplex)*N);
		memset(pmom[i], 0, sizeof(doublecomplex)*N);
	}

	// set starting configuration
	startconf(pmat, ALPHATILDE, NUMMAT, MATRIX_SIZE, LIEGROUP, conf, &L);

	// parameters only important for 8MM
	C2=(LIEGROUP-1.0)/(2.0*LIEGROUP) * L*(L+LIEGROUP);
	C3=/*C2**/((LIEGROUP*1.0-2.0)*(2.0*L+LIEGROUP*1.0))/(2.0*LIEGROUP);

	// compute action for model
	S = action_function(pmat, ALPHATILDE, MATRIX_SIZE, NUMMAT, C2, C3, LAMBDA, &S0, &S1);

	printf("S=%f\tS0=%f\tS1=%f\n", S,S0,S1);

	// begin integration steps k
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

		// initialize momentum with random configuration; this is not changed anymore after that
		if(k==1)
		{
			genmom(pmom, NUMMAT, MATRIX_SIZE, mommult);
			P = compmom(pmom, NUMMAT, MATRIX_SIZE, mass);
			H = S + P;
			fprintf(out5, "%f\t%f\t%f\t%f\t%f\n",H, S, P, S0, S1);
		}

		/*
		 * choose integration routine: Leapfrog or Omelyan
		 */

		//Leapfrog

		//			deltamom_function(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE, mass);
		//			deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS, C2, C3, LAMBDA);
		//			deltamom_function(pmat, pmom, EPS/2.0, NUMMAT, MATRIX_SIZE, mass);

		//Omelyan

		deltamom_function(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, mass);
		deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS/2.0, C2, C3, LAMBDA);
		deltamom_function(pmat, pmom, EPS*(1-2*ksi), NUMMAT, MATRIX_SIZE, mass);
		deltaaction_function(pmat, pmom, NUMMAT, MATRIX_SIZE, ALPHATILDE, EPS/2.0, C2, C3, LAMBDA);
		deltamom_function(pmat, pmom, ksi*EPS, NUMMAT, MATRIX_SIZE, mass);

		// compute final value of Action, Momentum and Hamiltonian; save in file
		Send = action_function(pmat, ALPHATILDE, MATRIX_SIZE, NUMMAT, C2, C3, LAMBDA, &S0, &S1);
		Pend = compmom(pmom, NUMMAT, MATRIX_SIZE, mass);
		Hend = Send + Pend;
		fprintf(out5, "%f\t%f\t%f\t%f\t%f\n",Hend, Send, Pend, S0, S1);
	}
	printf("\n");
	m_startTime = time(NULL);
	localtime_r(&m_startTime, &tim);
	printf("SU(3) ground state energy = %lf\n",-((LIEGROUP-2) ? MATRIX_SIZE : N)*LIEGROUP*ALPHATILDE*ALPHATILDE*ALPHATILDE*ALPHATILDE*C2/12);
	printf("Finished at %d:%d:%d\t%d/%d\n", tim.tm_hour, tim.tm_min, tim.tm_sec, tim.tm_mday, tim.tm_mon+1);

	fclose(out5);

	for(i=0;i<NUMMAT;i++)
	{
		free(pmat[i]);
		free(pmom[i]);
	}

	return 0;
}
