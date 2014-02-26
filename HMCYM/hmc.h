/*
 * hmc.h
 *
 *  Created on: 22 Feb 2011
 *      Author: tkaltenbrunner
 */

#ifndef HMC_H_
#define HMC_H_

int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size);
int addmattrace(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size);
int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_randphicplxtrace(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_gaussmomcplx(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
int gen_gaussmomcplxtrace(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
double mom(doublecomplex *pmom[], int nummat, int size);
double momtrace(doublecomplex *pmom[], int nummat, int size);
double actionYM(doublecomplex *pmat[], int size, int nummat, double mass);
int deltaYM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double eps, double mass);
int Traces(doublecomplex **pmat, int size, int nummat, int loopnumber, char *pathout);
int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat);
int computeGdd(doublecomplex *out, doublecomplex **in, int nummat, int size);
#endif

