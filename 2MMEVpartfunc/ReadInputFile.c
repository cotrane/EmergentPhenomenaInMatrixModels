#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputFile(int *lnum, int *size, double *eps, double *lacc, double *uacc, int *therm, int *steps, int *step, char *pathout, int *EVstore,
		double *mass, double *g, double *mult, int *compComm, int *compC, int *compD, int *compTraces, int *compPhi, int *compAComm, int *offDiag,
		int *nummat)
{
	FILE *read;
	char config[1000];
	char filenameread[256];
	char *check;

	sprintf(filenameread, "../../2MM.txt");
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
		if(!strcmp(config, "NUMMAT\n"))
		{
			check = fgets(config, 1000, read);
			*nummat = atoi(config);
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
		if(!strcmp(config, "MASS\n"))
		{
			check = fgets(config, 1000, read);
			*mass = atof(config);
		}
		if(!strcmp(config, "ParG\n"))
		{
			check = fgets(config, 1000, read);
			*g = atof(config);
		}
		if(!strcmp(config, "MULT\n"))
		{
			check = fgets(config, 1000, read);
			*mult = atof(config);
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
		//include computation of EV Traces
		if(!strcmp(config, "compTraces\n"))
		{
			check = fgets(config, 1000, read);
			*compTraces= atoi(config);
		}
		//include computation of EV of Phi=X+iY and Modulus r = sqrt(X^2+Y^2)
		if(!strcmp(config, "compPhi\n"))
		{
			check = fgets(config, 1000, read);
			*compPhi= atoi(config);
		}
		if(!strcmp(config, "compAComm\n"))
		{
			check = fgets(config, 1000, read);
			*compAComm= atoi(config);
		}
		if(!strcmp(config, "compOffDiag\n"))
		{
			check = fgets(config, 1000, read);
			*offDiag = atoi(config);
		}
	}

	fclose(read);
	return 0;
}
