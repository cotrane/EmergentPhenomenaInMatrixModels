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
}structconst2;
int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, int numprocs, int rank);
int addmom(doublecomplex *pmom[], doublecomplex *pmat[], double eps, int nummat, int size, double kappa);
int gen_randphicplx(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int gen_gaussmomcplx(doublecomplex *pmom[], int NUMMAT, int MATRIX_SIZE, double a);
double action(doublecomplex *pmat[], int nummat, int size, double kappa);
double mom(doublecomplex *pmom[], int nummat, int size, double lambda);
double action3MM(doublecomplex *pmat[], double coupling, int size, int nummat);
double action8MM(doublecomplex *pmat[], double *trYM, double *trCS, double coupling, int size, int nummat, double C2, double C3, double lambda, double *S0, double *S1, int numprocs, int rank);
int deltaX3MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps);
int deltaX8MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps, double C2, double C3, double lambda, int numprocs, int rank);
int deltaS1(doublecomplex *pmat[], doublecomplex *tmpmat[], int nummat, int size, double eps, double C2, double C3, double lambda, double alpha, int numprocs, int rank);
double S1(doublecomplex *pmat[], int size, double C2, double C3, double lambda, double alpha, int nummat, int numprocs, int rank);
int findStructConst(structconst2 **StructConst, int lie);
int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat, int numprocs, int rank);
#endif /* HMC_H_ */
