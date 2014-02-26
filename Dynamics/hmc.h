/*
 * hmc.h
 *
 *  Created on: 22 Feb 2011
 *      Author: tkaltenbrunner
 */

#ifndef HMC_H_
#define HMC_H_
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, double mass);
int addmom1MM(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size);
int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_gaussmomcplx(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
double action1MM(doublecomplex *pmat[], int nummat, int size);
double mom(doublecomplex *pmom[], int nummat, int size, double lambda);
double action3MM(doublecomplex *pmat[], double coupling, int size, int nummat);
double action3MMwrap(doublecomplex *pmat[], double coupling, int size, int nummat, double C2, double C3, double LAMBDA,
		double *S0, double *S1);
double action8MM(doublecomplex *pmat[], double coupling, int size, int nummat, double C2, double C3, double lambda,
		double *S0, double *S1);
int time_rev(doublecomplex *pmom[], doublecomplex *pmat[], int size, double eps, double ksi, int steps, int nummat,
		double alphatilde, double C2, double C3, double lambda);
int deltaX3MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps);
int deltaX3MMwrap(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps,
		double C2, double C3, double LAMBDA);
int deltaX8MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps, double C2,
		double C3, double lambda);
int deltaS1(doublecomplex *pmat[], doublecomplex *tmpmat[], int nummat, int size, double eps, double C2, double C3,
		double lambda, double alpha);
double S1(doublecomplex *pmat[], int size, double C2, double C3, double lambda, double alpha, int nummat);
double Mass(doublecomplex *pmat[], int nummat, int size, double constant);
int deltaMass(doublecomplex *pout[], doublecomplex *pmat[], int nummat, int size, double constant);
double mom(doublecomplex *pmom[], int nummat, int size, double mass);
double action1MMwrap(doublecomplex *pmat[], double coupling, int size, int nummat, double C2, double C3, double LAMBDA,
		double *S0, double *S1);
int deltaX1MMwrap(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps,
		double C2, double C3, double kappa);
int addmat1MM(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, double mass);
int gen_gaussmomcplxtrace(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
int gen_randphicplxtrace(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
#endif /* HMC_H_ */
