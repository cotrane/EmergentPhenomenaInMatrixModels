#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>

/*
 * computes EV distribution from EV's stored in file; parameters are read in via ReadInputFile
 */

int ReadInputFile(char *pathin, char *pathout, int *MATRIX_SIZE, int *Position, int *steps, double *binwidth);

int MATRIX_SIZE, position, steps;
double binwidth;

int main()
{
	int i,j,count=0, countlines=0;
	double RangeEV=0, interval=0, LargEV=0, SmallEV=0;
	double *eigenv;
	int *pev, Bin, BinsEV=0;
	char stream[30], pathin[1000], pathout[1000], *check;
	FILE *in, *out;
	char *filenamein, *filenameout;

	ReadInputFile(pathin, pathout, &MATRIX_SIZE, &position, &steps, &binwidth);
	if(position==0)
		position=1;

	filenamein=(char*) malloc(sizeof(char)*256);
	filenameout=(char*) malloc(sizeof(char)*256);
	eigenv = (double*) malloc(sizeof(double)*steps*MATRIX_SIZE);
	memset(eigenv, 0, sizeof(double)*steps*MATRIX_SIZE);

	sprintf(filenamein, "%s",pathin);
	sprintf(filenameout, "%s", pathout);
	in = fopen(filenamein,"r");
	if(in == 0){
		printf("\n==============================================================================\n");
		printf("wrong file name! %s\n", pathin);
		printf("==============================================================================\n");
	}
	out = fopen(filenameout, "w");

	count=0;
	countlines=1;
	check=fgets(stream, 30, in);
	while(check!=NULL)
	{
		if(countlines==(position*MATRIX_SIZE))
		{
			while(count!=((steps*MATRIX_SIZE)))
			{
				if(check==NULL){
					printf("just %d steps saved in file!!\nDistribution calculated with %dsteps\n", (count+position+MATRIX_SIZE)/MATRIX_SIZE, (count+MATRIX_SIZE)/(MATRIX_SIZE));
					break;
				}
				eigenv[count] = atof(stream);
				if(eigenv[count]>LargEV)
				{
					LargEV=eigenv[count];
				}
				else if(eigenv[count]<SmallEV)
				{
					SmallEV=eigenv[count];
				}
				if(count==0)
					printf("%d\t%f\n", count,eigenv[count]);
				count++;
				check=fgets(stream, 30, in);
			}
			break;
		}
		check=fgets(stream, 30, in);
		countlines++;
	}

	RangeEV = LargEV - SmallEV;
	BinsEV = RangeEV/binwidth +1;
	pev = (int*) malloc(BinsEV * sizeof(int));
	memset(pev, 0, BinsEV * sizeof(int));

	for(i=0;i<count;i++)
	{
		Bin = (eigenv[i]-SmallEV)/(binwidth*1.0);
		if(Bin > BinsEV)
		{
			printf("somethings wrong!\n");
		}
		if(Bin==BinsEV)
			Bin = BinsEV-1;
		pev[Bin] += 1;
	}

	for(j=0; j<BinsEV;j++)
	{
		interval = SmallEV + binwidth*j;
		fprintf(out, "%f\t%f\n", interval, pev[j]/(count*binwidth*1.0));
		fflush(out);
	}

	free(filenamein);
	free(filenameout);
	fclose(out);
	fclose(in);

	return 0;
}
