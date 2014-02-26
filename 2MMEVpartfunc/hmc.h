/*
 * hmc.h
 *
 *  Created on: 22 Feb 2011
 *      Author: tkaltenbrunner
 */

#ifndef HMC_H_
#define HMC_H_

#include <f2c.h>

typedef struct{
	int a;
	int b;
	int c;
	double val;
	int numEntries;
}structconst;
int addmat(double pmat[], double pmom[], double eps, int size);
int addmattrace(double pmat[], double pmom[], double eps, int size);
int gen_randphicplx(double pmat[], int MATRIX_SIZE);
int gen_randphicplxtrace(double pmat[], int MATRIX_SIZE, double mult);
int gen_gaussmomcplxtrace(double pmom[], int MATRIX_SIZE, double a);
int gen_gaussmomcplx(double pmom[], int MATRIX_SIZE, double a);
int gen_gaussMat_cplx_trace(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
double mom(double pmom[], int size);
double momtrace(double pmom[], int size);
double actionYM(double pmat[], int size, double mass);
int deltaYM(double pmat[], double pmom[], int size, double eps, double mass);
double actionYM_N(double pmat[], int size, double g);
int deltaYM_N(double pmat[], double pmom[], int size, double eps, double g);
int deltaYMtraceless(double pmat[], double pmom[], int size, double eps, double mass);
int GenGaussMatrixY(doublecomplex *pout, double *pmat, double mass, int size);
int ConstructAdjointAction(doublecomplex **out, doublecomplex **pmat, int size, int nummat);
int findStructConst(structconst **StructConst, int lie);
int Traces(doublecomplex **pmat, int size, int nummat, int loopnumber, char *pathout);
int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat);
#endif

