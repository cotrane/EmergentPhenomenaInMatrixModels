#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<string.h>

/*
 * reads in parameter supplied in command line via python script
 */
int ReadInputCL(int argc, char *argv[], int *step, int *size, double *alpha, int *loops, int *pos, char *pathout, char *pathin)
{
	int argNum=0;

	for(argNum = 1; argNum < argc; argNum++)
	{
		if(!strcmp(argv[argNum], "-par"))
		{
			if(!strcmp(argv[argNum + 1], "step"))
			{
				*step = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			else if(!strcmp(argv[argNum + 1], "loops"))
			{
				*loops = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			else if(!strcmp(argv[argNum + 1], "Position"))
			{
				*pos = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			else if(!strcmp(argv[argNum + 1], "MATRIX SIZE"))
			{
				*size = atoi(argv[argNum + 2]);
				argNum += 2;
			}
			else if(!strcmp(argv[argNum + 1], "ALPHA"))
			{
				*alpha = atof(argv[argNum + 2]);
				argNum += 2;
			}
			else if(!strcmp(argv[argNum + 1], "PATHIN"))
			{
				strcpy(pathin, argv[argNum + 2]);
				argNum += 2;
			}
			else if(!strcmp(argv[argNum + 1], "PATHOUT"))
			{
				strcpy(pathout, argv[argNum + 2]);
				argNum += 2;
			}
		}
	}

	return 0;
}
