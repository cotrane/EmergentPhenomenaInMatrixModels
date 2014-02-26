

#ifndef MONTECARLO_H_
#define MONTECARLO_H_
int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step,
		char *pathout, int *EVstore, int *nummat, int *liegroup, double *mass, int *compC);
int ReadInputFile(int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, char *pathout, int *EVstore, int *nummat
		, int *liegroup, double *mass, int *compC);
#endif /* MONTECARLO_H_ */
