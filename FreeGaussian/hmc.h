/*
 * hmc.h
 *
 *  Created on: 22 Feb 2011
 *      Author: tkaltenbrunner
 */

#ifndef HMC_H_
#define HMC_H_
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size);
int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_randphicplxtrace(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_gaussmomcplx(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
int gen_gaussmomcplxtrace(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
int gen_randphicplxdiagtrace(doublecomplex pmat[], int MATRIX_SIZE, double mult);
int GenGaussMatrix(doublecomplex *pout, double *pmat, double mass, int size);
int addmom(doublecomplex *pmom[], doublecomplex *pmat[], double eps, int nummat, int size, double mass);
double action(doublecomplex *pmat[], int nummat, int size, double mass);
double mom(doublecomplex *pmom[], int nummat, int size);
int genNonHermitianMatrix(doublecomplex *pmat, int size);
#endif

