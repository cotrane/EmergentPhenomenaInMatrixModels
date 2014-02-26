/*
 * MatrixMan.h
 *
 *  Created on: 26 Jul 2011
 *      Author: tkaltenbrunner
 */


#ifndef MATRIXMAN_H_
#define MATRIXMAN_H_
#include<f2c.h>

int printmat(char *use, doublecomplex *mat, int size);
int printmat_transp(char *use, doublecomplex *mat, int size);
int printmat_mathematica(char *use, doublecomplex *mat, int size);
int printmatdiff(doublecomplex *mat, doublecomplex *mat2, int size);
int AddMat(doublecomplex *out, doublecomplex *pcomm, doublecomplex *pmat, int size, double con1, double con2);
int AddMat2(doublecomplex *out, doublecomplex *pcomm, doublecomplex *pmat, int size, double con1, double con2);
int AddMat4(doublecomplex *out, doublecomplex *pcomm, int size, double con);
int AddMat4Dyn(doublecomplex *out, doublecomplex *pcomm, int size, double con);
int AddMat6(doublecomplex *out, doublecomplex *mat1, double con, int size);
int AddMat7(doublecomplex *out, doublecomplex *mat1, double con, int size);
int AddMat8(doublecomplex *out, doublecomplex *mat1, double con, int size);
int AddMat9(doublecomplex *out, doublecomplex *pmat, doublecomplex *pmat2, int size, double con);
void AntiComm(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double con);
void AntiComm3(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double con);
void AntiComm4(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double con);
int Comm(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con);
int Comm2(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con);
int Comm3(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con);
int Commend2(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con);
int diagMulti(double *out, doublecomplex *mat1, doublecomplex *mat2, int size, double constant);
int Multi5(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double constant);
int Multi6(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double constant);
void Multilambda(doublecomplex *out, doublecomplex *mat, doublecomplex *lambda, int size, int sizelambda);
int diagMultiCplx(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double constant);
int diagMultiI(double *out, doublecomplex *mat1, doublecomplex *mat2, int size, double constant);
int Inverse_posDef(doublecomplex *out, doublecomplex *in, int size);
int Inverse_indef(doublecomplex *out, doublecomplex *in, integer size);
#endif /* MATRIXMAN_H_ */
