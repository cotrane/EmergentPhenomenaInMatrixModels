#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include <time.h>

/*
 * computes mat[a]*mat[b] + mat[b]*mat[a]; anticommutator -> result is hermitian; saves computation time; does not include traceless condition!
 * used in GammaMatrices, 2MMEVpartfunc, Dynamics(8MM)
 */
void AntiComm(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double con)
{
	int i,j,k;

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			for(k=0; k<size; k++)
			{
				out[i*size+j].r += con*(mat[i*size+k].r * mat2[k*size+j].r - mat[i*size+k].i * mat2[k*size+j].i +
						mat2[i*size+k].r * mat[k*size+j].r - mat2[i*size+k].i * mat[k*size+j].i);
				out[i*size+j].i += con*(mat[i*size+k].r * mat2[k*size+j].i + mat[i*size+k].i * mat2[k*size+j].r +
						mat2[i*size+k].r * mat[k*size+j].i + mat2[i*size+k].i * mat[k*size+j].r);
			}
		}
	}
}

/*
 * used in Dynamics(8MM), HMC(8MM)
 */
void AntiComm3(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double con)
{
	int i,j,k;

	for(i=0; i<size; i++)
	{
		for(j=(i+1); j<size; j++)
		{
			for(k=0; k<size; k++)
			{
				out[i*size+j].r += con*(mat[j*size+k].r * mat2[k*size+i].r - mat[j*size+k].i * mat2[k*size+i].i +
						mat2[j*size+k].r * mat[k*size+i].r - mat2[j*size+k].i * mat[k*size+i].i);
				out[i*size+j].i += con*(mat[j*size+k].r * mat2[k*size+i].i + mat[j*size+k].i * mat2[k*size+i].r +
						mat2[j*size+k].r * mat[k*size+i].i + mat2[j*size+k].i * mat[k*size+i].r);
			}
			out[j*size+i].r = out[i*size+j].r;
			out[j*size+i].i = -out[i*size+j].i;
		}
		for(k=0;k<size;k++)
		{
			out[i*size+i].r += con*(mat[i*size+k].r * mat2[k*size+i].r - mat[i*size+k].i * mat2[k*size+i].i +
					mat2[i*size+k].r * mat[k*size+i].r - mat2[i*size+k].i * mat[k*size+i].i);
			out[i*size+i].i += con*(mat[i*size+k].r * mat2[k*size+i].i + mat[i*size+k].i * mat2[k*size+i].r +
					mat2[i*size+k].r * mat[k*size+i].i + mat2[i*size+k].i * mat[k*size+i].r);
		}
	}
}

// used in HMC(8MM),Dynamics(8MM)
void AntiComm4(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double con)
{
	int i,j,k;

	for(i=0; i<size; i++)
	{
		for(j=(i+1); j<size; j++)
		{
			for(k=0; k<size; k++)
			{
				out[i*size+j].r += con*(mat[i*size+k].r * mat2[k*size+j].r - mat[i*size+k].i * mat2[k*size+j].i + mat2[i*size+k].r * mat[k*size+j].r - mat2[i*size+k].i * mat[k*size+j].i);
				out[i*size+j].i += con*(mat[i*size+k].r * mat2[k*size+j].i + mat[i*size+k].i * mat2[k*size+j].r + mat2[i*size+k].r * mat[k*size+j].i + mat2[i*size+k].i * mat[k*size+j].r);
			}
			out[j*size+i].r = out[i*size+j].r;
			out[j*size+i].i = -out[i*size+j].i;
		}
		for(k=0;k<size;k++)
		{
			out[i*size+i].r += con*(mat[i*size+k].r * mat2[k*size+i].r - mat[i*size+k].i * mat2[k*size+i].i + mat2[i*size+k].r * mat[k*size+i].r - mat2[i*size+k].i * mat[k*size+i].i);
			out[i*size+i].i += con*(mat[i*size+k].r * mat2[k*size+i].i + mat[i*size+k].i * mat2[k*size+i].r + mat2[i*size+k].r * mat[k*size+i].i + mat2[i*size+k].i * mat[k*size+i].r);
		}
	}
}
