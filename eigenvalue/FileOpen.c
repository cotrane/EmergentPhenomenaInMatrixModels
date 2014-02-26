#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

/*
 * generates file in which EV's or EV distribution is saved and pointer to it; Save*-files use this pointer to copy the distribution into the file
 *
 * if operation==0: the eigenvalue distribution is stored
 * if operation==1: the individual eigenvalues are saved
 */

char *pfilename;
FILE *pout;
char *pfilenamest;
FILE *poutst;

void FileOpen(int num, int MATRIX_SIZE, int LOOPNUMBER, int operation, char *pathout, char *name, int parsnum, char *parsname[], double *parsval)
{
	int j;
	char tmp[256], filestring[1000];

	pout = (FILE*) malloc(sizeof(FILE));
	poutst = (FILE*) malloc(sizeof(FILE));

	if(operation==0)
	{
		pfilename = (char*) malloc(sizeof(char)*256);
	}
	if(operation==1)
	{
		pfilenamest = (char*) malloc(sizeof(char)*256);
	}

	memset(filestring, 0, sizeof(char)*1000);
	for(j=0;j<parsnum;j++)
	{
		if(!parsval[j])
			continue;
		else
		{
			sprintf(tmp, "-%s=%.2f", parsname[j], parsval[j]);
			strcat(filestring, tmp);
		}
	}

	if(operation == 0)
	{
			sprintf(pfilename, "%s%s%d-%dx%d%s.txt", pathout, name, num+1, MATRIX_SIZE, LOOPNUMBER, filestring);
			pout = fopen(pfilename,"w");
	}
	if(operation == 1)
	{
			sprintf(pfilenamest, "%s%s%dstore-%dx%d%s.txt", pathout, name, num+1, MATRIX_SIZE, LOOPNUMBER, filestring);
			poutst = fopen(pfilenamest,"a");
	}
}
