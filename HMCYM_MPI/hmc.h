/*
 * hmc.h
 *
 *  Created on: 22 Feb 2011
 *      Author: tkaltenbrunner
 */

#ifndef HMC_H_
#define HMC_H_
typedef struct{
	int a;
	int b;
	int c;
	double val;
	int numEntries;
}structconst;
double gauss_rand(int n);
double gauss_randvarmean(int n, double var, double mean);
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, int numprocs, int rank);
int addmattrace(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size);
int addmatalt(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size);
int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_randphicplxtrace(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_gaussmomcplx(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
int gen_gaussmomcplxtrace(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
double mom(doublecomplex *pmom[], int nummat, int size);
double momtrace(doublecomplex *pmom[], int nummat, int size);
double actionYM(doublecomplex *pmat[], int size, int nummat, double mass, int numprocs, int rank);
int deltaYM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double eps, double mass, int numprocs, int rank);
double actionYM_N(double pmat[], int size, double g);
int deltaYM_N(double pmat[], double pmom[], int size, double eps, double g);
int GenGaussMatrix(doublecomplex *pout, double *pmat, double mass, int size);
int findStructConst(structconst **StructConst, int lie);
int computeGdd(doublecomplex *out, doublecomplex **in, int nummat, int size, int rank, int numprocs);
int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat, int numprocs, int rank);
#endif

