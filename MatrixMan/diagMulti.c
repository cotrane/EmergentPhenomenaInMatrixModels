#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include <time.h>

/*
 * compute the trace of matrix multiplication
 * used in 2MM, Dynamics(3MM,8MM),RepsSUn,HMC,FreeGaussian,HMCYM
 */
int diagMulti(double *out, doublecomplex *mat1, doublecomplex *mat2, int size, double constant)
{
	int i,t;

	for(i = 0; i < size ;i++)
	{
		for(t=0; t < size; t++)
		{
			*out += constant*(mat1[i*size + t].r * mat2[t*size + i].r - mat1[i*size + t].i * mat2[t*size + i].i);
			*out += constant*(mat1[i*size + t].r * mat2[t*size + i].i + mat1[i*size + t].i * mat2[t*size + i].r);
		}
	}

	return 0;
}

/*
 * compute trace of matrix multiplication of non-Hermitian matrix; trace will be complex
 * used in HMCYM
 */
int diagMultiCplx(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double constant)
{
	int i,t;

	for(i = 0; i < size ;i++)
	{
		for(t=0; t < size; t++)
		{
			(*out).r += constant*(mat1[i*size + t].r * mat2[t*size + i].r - mat1[i*size + t].i * mat2[t*size + i].i);
			(*out).i += constant*(mat1[i*size + t].r * mat2[t*size + i].i + mat1[i*size + t].i * mat2[t*size + i].r);
		}
	}

	return 0;
}


/*
 * computes trace of 2 matrices where the 2nd matrix has an "i" as a coefficient!!!
 */
int diagMultiI(double *out, doublecomplex *mat1, doublecomplex *mat2, int size, double constant)
{
	int i,t;

	for(i = 0; i < size ;i++)												//multiplying matrices
	{
		for(t=0; t < size; t++)
		{
			*out += constant*(-mat1[i*size + t].r * mat2[t*size + i].i - mat1[i*size + t].i * mat2[t*size + i].r);
			*out += constant*(mat1[i*size + t].r * mat2[t*size + i].r - mat1[i*size + t].i * mat2[t*size + i].i);
		}
	}

	return 0;
}
