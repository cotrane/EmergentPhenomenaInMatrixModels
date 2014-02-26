#include<f2c.h>
#include<stdio.h>



/*
 * careful!! result of multiplication is NOT hermitian!! don't try to simplify it!!! :P
 * used in 2MM, Dynamics(8MM), HMC(8MM)
 */
int Multi5(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double constant)
{
	int i,j,k;

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			for(k=0; k<size; k++)
			{
				out[i*size+j].r += constant * (mat[i*size+k].r * mat2[k*size+j].r - mat[i*size+k].i * mat2[k*size+j].i);
				out[i*size+j].i += constant * (mat[i*size+k].r * mat2[k*size+j].i + mat[i*size+k].i * mat2[k*size+j].r);
			}
		}
	}

	return 0;
}

/*
 * same as Multi5 but adds (j,i)'th element to (i,j)
 * used in Dynamics(8MM), HMC(8MM)
 */
int Multi6(doublecomplex *out, doublecomplex *mat, doublecomplex *mat2, int size, double constant)
{
	int i,j,k;

	for(i=0; i<size; i++)
	{
		for(j=0; j<size; j++)
		{
			for(k=0; k<size; k++)
			{
				out[i*size+j].r += constant * (mat[j*size+k].r * mat2[k*size+i].r - mat[j*size+k].i * mat2[k*size+i].i);
				out[i*size+j].i += constant * (mat[j*size+k].r * mat2[k*size+i].i + mat[j*size+k].i * mat2[k*size+i].r);
			}
		}
	}

	return 0;
}
