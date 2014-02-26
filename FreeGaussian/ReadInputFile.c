#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputFile(int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, char *pathout, int *EVstore, int *nummat
		, int *liegroup, double *mass, int *compC)
{
	FILE *read;
	char config[1000];
	char filenameread[256];
	char *check;

	sprintf(filenameread, "../../FreeGaussian.txt");
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
		if(!strcmp(config, "PATHOUT\n"))
		{
			check = fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathout, config);
		}
		if(!strcmp(config, "EVstore\n"))
		{
			check = fgets(config, 1000, read);
			*EVstore = atoi(config);
		}
		if(!strcmp(config, "NUMMAT\n"))
		{
			check = fgets(config, 1000, read);
			*nummat = atoi(config);
		}
		if(!strcmp(config, "LIEGROUP\n"))
		{
			check = fgets(config, 1000, read);
			*liegroup = atoi(config);
		}
		if(!strcmp(config, "MASS\n"))
		{
			check = fgets(config, 1000, read);
			*mass = atof(config);
		}
		if(!strcmp(config, "compC\n"))
		{
			check = fgets(config, 1000, read);
			*compC = atoi(config);
		}
	}

	fclose(read);
	return 0;
}
