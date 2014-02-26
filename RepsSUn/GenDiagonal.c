#include <stdio.h>
#include <stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>

/*
 * generates the diagonal elements of the lie generators (element of the Cartan-subalgebra)
 */

int GenDiagonal(doublecomplex **ppmat, double **pplambda, int **ppstate, int **ppMultiplet, int **ppMultinum, int N, int dim, double a)
{
	int i,j,k,l=2,row,col,sum1,num;
	double entry, *norm;
	doublecomplex *pcart;

	norm = (double*) malloc(sizeof(double)*(N-1));

	//generating diagonal lambda's of fundamental rep to get diagonal operators for SU(N)
	for(i=2;i<(N+1);i++)
	{
		sum1=0;
		for(j=0;j<i;j++)
		{
			if(j<(i-1))
			{
				*(pplambda[i-2]+j) = 1.0;
				sum1+=1.0;
			}
			else
			{
				*(pplambda[i-2]+j) = -sum1;
			}
		}
	}
	//calculating their normalization
	for(i=0;i<(N-1);i++)
	{
		norm[i]=0;
		for(j=0;j<N;j++)
		{
			norm[i] += *(pplambda[i]+j)* *(pplambda[i]+j);
		}
		norm[i] /= 2.0;
		norm[i] = 1.0/sqrt(norm[i]);
	}

	for(i=((N*N-1)-(N-1));i<(N*N-1);i++)
	{
		num=l*l-2;
		pcart = ppmat[num];
		for(j=0;j<dim;j++)
		{
			entry=0;
			for(k=0;k<N;k++)
			{
				entry+= a/2.0* *(pplambda[i-((N*N-1)-(N-1))]+k)*norm[i-((N*N-1)-(N-1))]* *(ppstate[j]+k);
			}
			row=abs(j-(dim-1));
			col=abs(j-(dim-1));
			pcart[j*dim+j].r = entry;
		}
		l++;
	}



	return 0;
}
