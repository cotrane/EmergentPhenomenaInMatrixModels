/*
 * hmc.h
 *
 *  Created on: 22 Feb 2011
 *      Author: tkaltenbrunner
 */

#ifndef HMC_H_
#define HMC_H_

int gen_randphicplx(doublecomplex *pmat[], int MATRIX_SIZE);
int gen_gaussmomcplx(doublecomplex *pmom[], int MATRIX_SIZE, double ac);
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, double lambda);
int addmom(doublecomplex *pmom[], doublecomplex *pmat[], double eps, int nummat, int size, double kappa);
double action(doublecomplex *pmat[], int nummat, int size, double kappa);
double mom(doublecomplex *pmom[], int nummat, int size, double lambda);
double action3MM(doublecomplex *pmat[], double *trYM1, double *trCS1, double coupling, int size, int nummat);
double action3MMwrap(doublecomplex *pmat[], double *trYM1, double *trCS1, double coupling, int size, int nummat,
		double C2, double C3, double LAMBDA, double *S0, double *S1);
double action8MM(doublecomplex *pmat[], double *trYM1, double *trCS1, double coupling, int size, int nummat,
		double C2, double C3, double lambda, double *S0, double *S1);
int deltaX3MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps);
int deltaX3MMwrap(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps,
		double C2, double C3, double LAMBDA);
int deltaX8MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps,
		double C2, double C3, double lambda);
int deltaS1(doublecomplex *pmat[], doublecomplex *tmpmat[], int nummat, int size, double eps, double C2, double C3,
		double lambda, double alpha);
double S1(doublecomplex *pmat[], int size, double C2, double C3, double lambda, double alpha, int nummat);
int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat);
typedef struct{
	int a;
	int b;
	int c;
	double val;
	int numEntries;
}structconst2;
int findStructConst(structconst2 **StructConst, int lie);
#endif

