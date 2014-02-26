#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include <time.h>

/*
 * Comm where (i,j) is added to (j,i); using hermiticity and including tracelessness cond.
 * used in Dynamics(8MM), HMCYM, HMC(8MM), HMC(3MM),HMCYM
 */
int Commend2(doublecomplex *out, doublecomplex *mat1, doublecomplex *mat2, int size, double con)
{
	int i,j,t;
	double lastel=0;

	for(i = 0; i < size ;i++)												//calculating commutators
	{
		for(j=(i+1); j < size; j++)
		{
			for(t=0; t < size; t++)
			{
				out[i*size + j].r += con*(mat1[j*size + t].r * mat2[t*size + i].r - mat1[j*size + t].i * mat2[t*size + i].i - mat2[j*size + t].r * mat1[t*size + i].r + mat2[j*size + t].i * mat1[t*size + i].i);
				out[i*size + j].i += con*(mat1[j*size + t].r * mat2[t*size + i].i + mat1[j*size + t].i * mat2[t*size + i].r - mat2[j*size + t].r * mat1[t*size + i].i - mat2[j*size + t].i * mat1[t*size + i].r);
				if(j==(i+1))
				{
					out[i*size + i].r += con*(mat1[i*size + t].r * mat2[t*size + i].r - mat1[i*size + t].i * mat2[t*size + i].i - mat2[i*size + t].r * mat1[t*size + i].r + mat2[i*size + t].i * mat1[t*size + i].i);
					lastel += con*(mat1[i*size + t].r * mat2[t*size + i].r - mat1[i*size + t].i * mat2[t*size + i].i - mat2[i*size + t].r * mat1[t*size + i].r + mat2[i*size + t].i * mat1[t*size + i].i);
				}
			}
			out[j*size + i].r = out[i*size + j].r;
			out[j*size + i].i = -out[i*size + j].i;
		}
	}
	out[(size-1)*size+(size-1)].r += -lastel;

	return 0;
}
