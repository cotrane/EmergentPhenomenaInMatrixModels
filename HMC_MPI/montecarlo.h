/*
 * montecarlo.h
 *
 *  Created on: 05.07.2010
 *      Author: Thomas
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include<stdio.h>

int thecode(int mat, double *mat_c);
int symmsc(double *mat_c2);
int symmsc2(double *mat_c2, int mat);
int structconst(double *mat_c);
int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, int *steps,
		int *step, double *lambda, char *pathout, int *config, int *EVstore, double *mommult, int *numrep, int **rep, int *diag, double *prop1,
		double *prop2, int *compC, int *compD, int *compEV, int *compComm, int *compComm2, int *compB);
int ReadInputFile(int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, double *lambda,
		char *pathout, int *conf, double *mommult, int *EVstore, int *numrep, int **rep, int *diag, double *prop1, double *prop2, int *compC,
		int *compD, int *compEV, int *compComm, int *compComm2, int *compB);
int startconf(doublecomplex *pmat[], double alpha, int nummat, int size, int liegroup, int numrep, int *rep, int conf, int *L, int diag, double prop1, double prop2);
#endif /* MONTECARLO_H_ */
