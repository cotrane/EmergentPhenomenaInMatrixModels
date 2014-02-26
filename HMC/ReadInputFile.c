#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>
#include "model.h"

/*
 * reading in parameters from file HMC.txt
 */

int ReadInputFile(int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step,
		double *lambda, char *pathout, int *conf, double *mommult, int *EVstore, int *numrep, int **rep, int *diag, double *prop1, double *prop2,
		int *compEV, int *compComm, int *compC, int *compD, int *compB)
{
	int i;
	FILE *read;
	char config[1000], *repchar=0;
	char filenameread[256];
	char *check;

	sprintf(filenameread, "../../HMC.txt");
	read = fopen(filenameread, "r");
	while(fgets(config, 1000, read)!= NULL)
	{
		if(!strcmp(config, "LOOPNUMBER\n"))
		{
			check = fgets(config, 1000, read);
			*lnum = atoi(config);
		}
		if(!strcmp(config, "SIZE\n"))
		{
			check = fgets(config, 1000, read);
			*size = atoi(config);
		}
		if(!strcmp(config, "ALPHATILDE\n"))
		{
			check = fgets(config, 1000, read);
			*atilde = atof(config);
		}
		if(!strcmp(config, "EPS\n"))
		{
			check = fgets(config, 1000, read);
			*eps = atof(config);
		}
		if(!strcmp(config, "LACC\n"))
		{
			check = fgets(config, 1000, read);
			*lacc = atof(config);
		}
		if(!strcmp(config, "UACC\n"))
		{
			check = fgets(config, 1000, read);
			*uacc = atof(config);
		}
		if(!strcmp(config, "THERM\n"))
		{
			check = fgets(config, 1000, read);
			*therm = atoi(config);
		}
		if(!strcmp(config, "STEPS\n"))
		{
			check = fgets(config, 1000, read);
			*steps = atoi(config);
		}
		if(!strcmp(config, "STEP\n"))
		{
			check = fgets(config, 1000, read);
			*step = atoi(config);
		}
		if(!strcmp(config, "LAMBDA\n"))
		{
			check = fgets(config, 1000, read);
			*lambda = atof(config);
		}
		if(!strcmp(config, "PATHOUT\n"))
		{
			check = fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathout, config);
		}
		if(!strcmp(config, "START\n"))
		{
			check = fgets(config, 1000, read);
			*conf = atoi(config);
		}
		if(!strcmp(config, "MOMMULT\n"))
		{
			check = fgets(config, 1000, read);
			*mommult = atof(config);
		}
		if(!strcmp(config, "EVstore\n"))
		{
			check = fgets(config, 1000, read);
			*EVstore = atoi(config);
		}
		if(!strcmp(config, "REPS\n"))
		{
			check = fgets(config, 1000, read);
			int curr_char = 0;
			while(config[curr_char] != '\n')
			{
				if(config[curr_char] == '-')
					(*numrep)++;
				curr_char++;
			}
			(*numrep)++;
			*rep = (int*) malloc(sizeof(int)*(*numrep));
			memset(*rep, 0, sizeof(int)*(*numrep));
			repchar = strtok(config, "-");
			(*rep)[0] = atoi(repchar);
			for(i=1;i<(*numrep);i++)
			{
				repchar = 0;
				repchar = strtok(NULL, "-");
				(*rep)[i] = atoi(repchar);
			}
		}
		if(!strcmp(config, "DIAG\n"))
		{
			check = fgets(config, 1000, read);
			*diag = atoi(config);
		}
		if(!strcmp(config, "PROP1\n"))
		{
			check = fgets(config, 1000, read);
			*prop1 = atof(config);
		}
		if(!strcmp(config, "PROP2\n"))
		{
			check = fgets(config, 1000, read);
			*prop2 = atof(config);
		}
		if(!strcmp(config, "compEV\n"))
		{
			check = fgets(config, 1000, read);
			*compEV = atoi(config);
		}
		if(!strcmp(config, "compComm\n"))
		{
			check = fgets(config, 1000, read);
			*compComm = atoi(config);
		}
		if(!strcmp(config, "compC\n"))
		{
			check = fgets(config, 1000, read);
			*compC = atoi(config);
		}
		if(!strcmp(config, "compD\n"))
		{
			check = fgets(config, 1000, read);
			*compD = atoi(config);
		}
		if(!strcmp(config, "compB\n"))
		{
			check = fgets(config, 1000, read);
			*compB = atoi(config);
		}
	}

	fclose(read);
	return 0;
}
