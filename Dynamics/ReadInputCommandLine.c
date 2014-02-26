#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, double *lambda, char *pathout, int *conf, double *mommult, double *mass)
{
	int argNum=0;

	for(argNum = 1; argNum < argc; argNum++)
	{
		if(!strcmp(argv[argNum], "-par"))
		{
			if(!strcmp(argv[argNum + 1], "LOOPNUMBER"))
			{
				*lnum = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "SIZE"))
			{
				*size = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "ALPHATILDE"))
			{
				*atilde = atof(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "EPS"))
			{
				*eps = atof(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "LACC"))
			{
				*lacc = atof(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "UACC"))
			{
				*uacc = atof(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "THERM"))
			{
				*therm = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "STEPS"))
			{
				*steps = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "STEP"))
			{
				*step = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "LAMBDA"))
			{
				*lambda = atof(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "PATHOUT"))
			{
				strcpy(pathout, argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "START"))
			{
				*conf = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "MOMMULT"))
			{
				*mommult = atof(argv[argNum + 2]);
				argNum += 2;
			}
			if(!strcmp(argv[argNum + 1], "MASS"))
			{
				*mass = atof(argv[argNum + 2]);
				argNum += 2;
			}
		}
	}

	return 0;
}
