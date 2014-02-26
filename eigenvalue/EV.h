/*
 * EV.h
 *
 *  Created on: 5 Nov 2010
 *      Author: tkaltenbrunner
 */

#ifndef EV_H_
#define EV_H_

int findEigenvalues(doublereal *pVec, doublecomplex *m_pMatrix, integer matrix_size);
int findEigenvalues_vectors(doublereal *pVec, doublecomplex *m_pMatrix, integer m_nDim);
void SaveEV_1vector(double pSavEV[], double *LargEV, double *SmallEV, doublecomplex *eigenvectors, int MATRIX_SIZE);
void SaveEV(double *pSavEV, double *LargEV, double *SmallEV, doublecomplex *pmat, int MATRIX_SIZE, int *numev);
void FileClose(int operation);
void FileOpen(int num, int MATRIX_SIZE, int LOOPNUMBER, int operation, char *pathout, char *name, int parsnum, char *parsname[], double *parsval);
void DataDistribution(double *pSavEV, double LargEV, double SmallEV, double BINWIDTH, int num, int numev);
void StoreData(double *pSavEV, int *numev);
int EVonthefly(int **out, doublecomplex *in, int size, double binwidth, int *bins, double *LargeEV,	double *SmallEV, int *numev);
int SaveEVendofflight(int *in, int BinsEV, double binwidth, double SmallEV, int numev);
int findEigenvalues_NonHermitian(doublecomplex *pVec, doublecomplex *m_pMatrix, integer m_nDim);
void SaveEV_nonHermitian(double *pSavEVr, double *pSavEVi, doublecomplex *LargEV, doublecomplex *SmallEV, doublecomplex *pmat,
		int MATRIX_SIZE, int *numev);
int EVtoDist_nonHermitian(int **out1, int **out2, doublecomplex *in, int size, double binwidth, int bins[2], doublecomplex *LargeEV,
		doublecomplex *SmallEV, int *numev);
int EVtoDist_Modulus_nH(int **out, doublecomplex *in, int size, double binwidth, int *bins,	double *LargeEVmod, double *SmallEVmod, int *numEVmod);

#endif /* EV_H_ */
