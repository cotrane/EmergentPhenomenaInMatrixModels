

#ifndef MONTECARLO_H_
#define MONTECARLO_H_
int ReadInputFile(int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, char *pathout, int *EVstore,
		double *mass, double *g, double *mult, int *compComm, int *compC, int *compD, int *compTraces, int *compPhi, int *compAComm, int *offDiag,
		int *nummat);
int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step,
		char *pathout, int *EVstore, double *mass, double *g, double *mult, int *compComm, int *compC, int *compD, int *compTraces, int *compPhi,
		int *compAComm, int *offDiag, int *nummat);
int printmat_transp(char *use, doublecomplex *mat, int size);
#endif /* MONTECARLO_H_ */
