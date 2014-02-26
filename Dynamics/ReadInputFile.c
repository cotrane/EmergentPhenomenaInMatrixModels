#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputFile(int *lnum, int *size, double *atilde, double *eps, double *lacc, double *uacc, int *therm,
		double *lambda, char *pathout, int *conf, double *mommult, double *mass)
{
	FILE *read;
	char config[1000];
	char filenameread[256];
	char *check;

	sprintf(filenameread, "../../Dynamics.txt");
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
		if(!strcmp(config, "MASS\n"))
		{
			check = fgets(config, 1000, read);
			*mass = atof(config);
		}
	}

	fclose(read);
	return 0;
}
