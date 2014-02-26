#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

extern char *pfilename;
extern FILE *pout;
extern char *pfilenamest;
extern FILE *poutst;

/*
 * closes file opened in FileOpen in which the EV distribution/ EV's are stored
 */
void FileClose(int operation)
{
	if(operation==0)
	{
		free(pfilename);
		fclose(pout);
	}
	if(operation==1)
	{
		free(pfilenamest);
		fclose(poutst);
	}
}
