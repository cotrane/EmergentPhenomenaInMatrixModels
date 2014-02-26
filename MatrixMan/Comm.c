#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include <time.h>

/*
 * output is anti-hermitian!!! signs!!!!
 * used in 2MM, Dynamics(3MM,8MM),HMCYM,FreeGaussian, HMC, HMCYM
 */
int Comm(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con)
{
	int i,j,t;

	for(i = 0; i < size ;i++)												//calculating commutators
	{
		for(j=(i+1); j < size; j++)
		{
			for(t=0; t < size; t++)
			{
				out[i*size + j].r += con*(mat1[i*size + t].r * mat2[t*size + j].r - mat1[i*size + t].i * mat2[t*size + j].i - mat2[i*size + t].r * mat1[t*size + j].r + mat2[i*size + t].i * mat1[t*size + j].i);
				out[i*size + j].i += con*(mat1[i*size + t].r * mat2[t*size + j].i + mat1[i*size + t].i * mat2[t*size + j].r - mat2[i*size + t].r * mat1[t*size + j].i - mat2[i*size + t].i * mat1[t*size + j].r);
			}
			out[j*size + i].r = -out[i*size + j].r;
			out[j*size + i].i = out[i*size + j].i;
		}
		for(t=0;t<size;t++)
		{
			out[i*size + i].i += con*(mat1[i*size + t].r * mat2[t*size + i].i + mat1[i*size + t].i * mat2[t*size + i].r - mat2[i*size + t].r * mat1[t*size + i].i - mat2[i*size + t].i * mat1[t*size + i].r);
		}
	}

	return 0;
}

/*
 * simple commutator with coefficient
 * used in Dynamics(8MM)
 */

int Comm2(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con)
{
	int i,j,t;

	for(i = 0; i < size ;i++)												//calculating commutators
	{
		for(j=0; j < size; j++)
		{
			for(t=0; t < size; t++)
			{
				out[i*size + j].r += con*(mat1[i*size + t].r * mat2[t*size + j].r - mat1[i*size + t].i * mat2[t*size + j].i - mat2[i*size + t].r * mat1[t*size + j].r + mat2[i*size + t].i * mat1[t*size + j].i);
				out[i*size + j].i += con*(mat1[i*size + t].r * mat2[t*size + j].i + mat1[i*size + t].i * mat2[t*size + j].r - mat2[i*size + t].r * mat1[t*size + j].i - mat2[i*size + t].i * mat1[t*size + j].r);
			}
		}
	}


	return 0;
}

/*
 * HERMITIAN OUTPUT
 * used in Dynamics(3MM), HMC(3MM)
 */
int Comm3(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con)
{
	int i,j,t;

	for(i = 0; i < size ;i++)												//calculating commutators
	{
		for(j=(i+1); j < size; j++)
		{
			for(t=0; t < size; t++)
			{
				out[i*size + j].r += con*(mat1[i*size + t].r * mat2[t*size + j].r - mat1[i*size + t].i * mat2[t*size + j].i - mat2[i*size + t].r * mat1[t*size + j].r + mat2[i*size + t].i * mat1[t*size + j].i);
				out[i*size + j].i += con*(mat1[i*size + t].r * mat2[t*size + j].i + mat1[i*size + t].i * mat2[t*size + j].r - mat2[i*size + t].r * mat1[t*size + j].i - mat2[i*size + t].i * mat1[t*size + j].r);
			}
			out[j*size + i].r = out[i*size + j].r;
			out[j*size + i].i = -out[i*size + j].i;
		}
		for(t=0;t<size;t++)
		{
			out[i*size + i].r += con*(mat1[i*size + t].r * mat2[t*size + i].r - mat1[i*size + t].i * mat2[t*size + i].i - mat2[i*size + t].r * mat1[t*size + i].r + mat2[i*size + t].i * mat1[t*size + i].i);
			out[i*size + i].i += con*(mat1[i*size + t].r * mat2[t*size + i].i + mat1[i*size + t].i * mat2[t*size + i].r - mat2[i*size + t].r * mat1[t*size + i].i - mat2[i*size + t].i * mat1[t*size + i].r);
		}
	}

	return 0;
}
