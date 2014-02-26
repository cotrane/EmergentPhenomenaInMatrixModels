#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step,
		char *pathout, int *EVstore, int *nummat, double *mass, int *start, int *compC, int *compD, int *compEV, int *compComm)
{
	int argNum=0;

	for(argNum = 1; argNum + 2 < argc; argNum++)
	{
		if(!strcmp(argv[argNum], "-par"))
		{
			if(!strcmp(argv[argNum + 1], "LOOPNUMBER"))
			{
				*lnum = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "SIZE"))
			{
				*size = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "EPS"))
			{
				*eps = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "LACC"))
			{
				*lacc = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "UACC"))
			{
				*uacc = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "THERM"))
			{
				*therm = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "STEPS"))
			{
				*steps = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "STEP"))
			{
				*step = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "PATHOUT"))
			{
				strcpy(pathout, argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "EVstore"))
			{
				*EVstore = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "NUMMAT"))
			{
				*nummat = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "MASS"))
			{
				*mass = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "START"))
			{
				*start = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "compC"))
			{
				*compC = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "compD"))
			{
				*compD = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "compEV"))
			{
				*compEV = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "compComm"))
			{
				*compComm = atoi(argv[argNum + 2]);
			}
			argNum += 2;
		}
	}

	return 0;
}
