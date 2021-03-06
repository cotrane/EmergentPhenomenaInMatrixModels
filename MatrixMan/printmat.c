#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<f2c.h>


/*
 * printf command for matrix in screen formatting plus it's trace
 */
int printmat(char *use, doublecomplex *mat, int size)
{
	int j,k;
	double trrand=0, trrandi=0;

	printf("\n");
	printf("%s\n", use);

	for(j=0;j<size;j++)
	{
		for(k=0;k<size;k++)
		{
			if(k<(size-1))
				printf("%lf %lf\t", mat[j*size+k].r /*- mat[k*size+j].r*/, mat[j*size+k].i /*+ mat[k*size+j].i*/);
			else
				printf("%lf %lf\n", mat[j*size+k].r  /*- mat[k*size+j].r*/, mat[j*size+k].i /*+ mat[k*size+j].i*/);
		}
	}
	printf("\n");


	for(j=0;j<size;j++)
	{
		trrand += mat[j*size+j].r;
		trrandi += mat[j*size+j].i;
	}
	printf("%.10lf\t%.10f\n", trrand, trrandi);

	return 0;
}

/*
 * printf command for matrix in mathematica-input format plus it's trace
 */
int printmat_mathematica(char *use, doublecomplex *mat, int size)
{
	int j,k;
	double trrand=0, trrandi=0;

	printf("\n");
	printf("%s\n", use);
																				//mathematica formatting
	printf("{");
	for(j=0;j<size;j++)
	{	printf("{");
		for(k=0;k<size;k++)
		{
			if(k<(size-1))
			{
				if(mat[j*size+k].r!=0)
				{
					printf("%.10lf", mat[j*size+k].r);
				}
				else if(mat[j*size+k].r==0)
				{
					printf("%.0lf", mat[j*size+k].r);
				}

				if(mat[j*size+k].i!=0)
				{
					printf("+(%lfI),", mat[j*size+k].i);
				}
				else if(mat[j*size+k].i==0)
				{
					printf(",");
				}
			}
			else
			{
				if(mat[j*size+k].r!=0)
				{
					printf("%lf", mat[j*size+k].r);
				}
				else if(mat[j*size+k].r==0)
				{
					printf("%.0lf", mat[j*size+k].r);
				}

				if(mat[j*size+k].i!=0)
				{
					printf("+(%lfI)},", mat[j*size+k].i);
				}
				else if(mat[j*size+k].i==0)
				{
					printf("},");
				}
			}
		}
	}
	printf("}\n");


	for(j=0;j<size;j++)
	{
		trrand += mat[j*size+j].r;
		trrandi += mat[j*size+j].i;
	}
	printf("%.10lf\t%.10f\n", trrand, trrandi);

	return 0;
}

/*
 * prints the transposed of matrix on screen in screen formatting plus it's trace
 */
int printmat_transp(char *use, doublecomplex *mat, int size)
{
	int j,k;
	double trrand=0, trrandi=0;

	printf("\n");
	printf("%s\n", use);

	for(j=0;j<size;j++)
	{
		for(k=0;k<size;k++)
		{
			if(k<(size-1))
				printf("%lf %lf\t", mat[k*size+j].r /*- mat[k*size+j].r*/, mat[k*size+j].i /*+ mat[k*size+j].i*/);
			else
				printf("%lf %lf\n", mat[k*size+j].r  /*- mat[k*size+j].r*/, mat[k*size+j].i /*+ mat[k*size+j].i*/);
		}
	}
	printf("\n");


	for(j=0;j<size;j++)
	{
		trrand += mat[j*size+j].r;
		trrandi += mat[j*size+j].i;
	}
	printf("%.10lf\t%.10f\n", trrand, trrandi);

	return 0;
}

/*
 * printf command for difference of two matrices
 */
int printmatdiff(doublecomplex *mat, doublecomplex *mat2, int size)
{
	int j,k;

																				//screen formatting
	for(j=0;j<size;j++)
	{
		for(k=0;k<size;k++)
		{
			if(k<(size-1))
				printf("%lf %lf\t", mat[j*size+k].r - mat2[j*size+k].r, mat[j*size+k].i - mat2[j*size+k].i);
			else
				printf("%lf %lf\n", mat[j*size+k].r  - mat2[j*size+k].r, mat[j*size+k].i - mat2[j*size+k].i);
		}
	}

	return 0;
}

