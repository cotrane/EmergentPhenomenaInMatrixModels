#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputCL(int argc, char *argv[], int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, int *steps,
		int *step, double *lambda, char *pathout, int *conf, int *EVstore, double *mommult, int *numrep, int **rep, int *diag, double *prop1,
		double *prop2, int *compC, int *compD, int *compEV, int *compComm, int *compComm2, int *compB)
{
	int argNum=0, i;
	char strrep[1000];
	char *repchar=0;

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
			if(!strcmp(argv[argNum + 1], "ALPHATILDE"))
			{
				*atilde = atof(argv[argNum + 2]);
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
			if(!strcmp(argv[argNum + 1], "LAMBDA"))
			{
				*lambda = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "PATHOUT"))
			{
				strcpy(pathout, argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "START"))
			{
				*conf = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "MOMMULT"))
			{
				*mommult = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "EVstore"))
			{
				*EVstore = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "REPS"))
			{
				strcpy(strrep, argv[argNum+2]);
				int curr_char = 0;
				while(strrep[curr_char] != 0)
				{
					if(strrep[curr_char] == '-')
						(*numrep)++;
					curr_char++;
				}
				(*numrep)++;
				*rep = (int*) malloc(sizeof(int)*(*numrep));
				memset(*rep, 0, sizeof(int)*(*numrep));
				repchar = strtok(argv[argNum+2], "-");
				(*rep)[0] = atoi(repchar);
				for(i=1;i<(*numrep);i++)
				{
					repchar = 0;
					repchar = strtok(NULL, "-");
					(*rep)[i] = atoi(repchar);
				}
			}
			if(!strcmp(argv[argNum + 1], "DIAG"))
			{
				*diag = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "PROP1"))
			{
				*prop1 = atof(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "PROP2"))
			{
				*prop2 = atof(argv[argNum + 2]);
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
			if(!strcmp(argv[argNum + 1], "compComm2"))
			{
				*compComm2 = atoi(argv[argNum + 2]);
			}
			if(!strcmp(argv[argNum + 1], "compB"))
			{
				*compB = atoi(argv[argNum + 2]);
			}
			argNum += 2;
		}
	}

	return 0;
}
