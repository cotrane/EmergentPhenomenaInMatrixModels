#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

/*
 * reads in parameters supplied via ErrorEval.txt
 */
int ReadInputFile(int *step, int *size, double *alpha, int *loops, int *position, char *pathout, char *pathin)
{
	char filenameread[1000], config[1000];
	FILE *read;

	sprintf(filenameread, "../../ErrorEval.txt");
	read = fopen(filenameread, "r");
	while(fgets(config, 1000, read)!= NULL)
	{
		if(!strcmp(config, "PATHIN\n"))
		{
			fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathin, config);
		}
		else if(!strcmp(config, "PATHOUT\n"))
		{
			fgets(config, 1000, read);
			config[strlen(config) - 1] = 0;
			strcpy(pathout, config);
		}
		else if(!strcmp(config, "MATRIX_SIZE\n"))
		{
			fgets(config, 1000, read);
			*size = atoi(config);
		}
		else if(!strcmp(config, "ALPHA\n"))
		{
			fgets(config, 1000, read);
			*alpha = atof(config);
		}
		else if(!strcmp(config, "Position\n"))
		{
			fgets(config, 1000, read);
			*position = atoi(config);
		}
		else if(!strcmp(config, "loops\n"))
		{
			fgets(config, 1000, read);
			*loops = atoi(config);
		}
		else if(!strcmp(config, "step\n"))
		{
			fgets(config, 1000, read);
			*step = atoi(config);
		}
	}

	return 0;
}
