#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

int ReadInputFile(char *pathin, char *pathout, int *MATRIX_SIZE, int *Position, int *steps, double *binwidth)
{
	FILE *read;
	char config[1000];
	char filenameread[256];
	char *check;

	sprintf(filenameread, "../../EVEvaluation.txt");
	read = fopen(filenameread, "r");
	while(fgets(config, 1000, read)!= NULL)
	{
		if(!strcmp(config, "PATHIN\n"))
		{
			check = fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathin, config);
		}
		if(!strcmp(config, "PATHOUT\n"))
		{
			check = fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathout, config);
		}
		if(!strcmp(config, "MATRIX_SIZE\n"))
		{
			check = fgets(config, 1000, read);
			*MATRIX_SIZE = atoi(config);
		}
		if(!strcmp(config, "Position\n"))
		{
			check = fgets(config, 1000, read);
			*Position = atoi(config);
		}
		if(!strcmp(config, "STEPS\n"))
		{
			check = fgets(config, 1000, read);
			*steps = atoi(config);
		}
		if(!strcmp(config, "BINWIDTH\n"))
		{
			check = fgets(config, 1000, read);
			*binwidth = atof(config);
		}
	}

	fclose(read);
	return 0;
}
