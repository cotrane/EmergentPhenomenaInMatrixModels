

#ifndef MONTECARLO_H_
#define MONTECARLO_H_
int printmat(char *use, doublecomplex *mat, int size);
int printmatdiff(doublecomplex *mat, doublecomplex *mat2, int size);
int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step,
		char *pathout, int *EVstore, int *nummat, double *mass, int *start, int *compC, int *compD, int *compEV, int *compComm,
		int *compTraces, int *compTraceX2, int *compGdd);
int ReadInputFile(int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, char *pathout, int *EVstore,
		int *nummat, double *mass, int *start, int *compC, int *compD, int *compEV, int *compComm, int *compTraces, int *compTraceX2,
		int *compGdd);
#endif /* MONTECARLO_H_ */
