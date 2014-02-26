/*
 * montecarlo.h
 *
 *  Created on: 05.07.2010
 *      Author: Thomas
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_
int thecode(int mat, double *mat_c);
int symmsc(double *mat_c2);
int symmsc2(double *mat_c2, int mat);
int structconst(double *mat_c);
int dtensor(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int computesc(doublecomplex *pmat[], int NUMMAT, int MATRIX_SIZE);
int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm,
		double *lambda, char *pathout, int *conf, double *mommult, double *mass);
int ReadInputFile(int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, double *lambda,
		char *pathout, int *conf, double *mommult, double *mass);
int startconf(doublecomplex *pmat[], double alpha, int nummat, int size, int liegroup, int conf, int *L);
#endif /* MONTECARLO_H_ */
