#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include<stdio.h>

/*
 * careful!!!!!! just for use with coefficient "i" in 2nd matrix !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * AddMat adds a doublecommutator with a commutator with im coefficient -> hermitian
 * used in Dynamics(3MM), HMC(3MM)
 */

int AddMat(doublecomplex *out, doublecomplex *pcomm, doublecomplex *pmat, int size, double con1, double con2)
{
	int j,k;

	for(j=0;j<size;j++)
	{
		for(k=(j+1);k<size;k++)
		{
			out[j*size + k].r += con1 * pcomm[k*size + j].r - con2 * pmat[k*size + j].i;
			out[j*size + k].i += con1 * pcomm[k*size + j].i + con2 * pmat[k*size + j].r;

			out[k*size + j].r = out[j*size + k].r;
			out[k*size + j].i = -out[j*size + k].i;
		}
		out[j*size + j].r += con1 * pcomm[j*size + j].r - con2 * pmat[j*size + j].i;
	}

	return 0;
}

/*
 * adding a commutator and a matrix -> antihermitian!!
 * used in HMC(3MM)
 */
int AddMat2(doublecomplex *out, doublecomplex *pcomm, doublecomplex *pmat, int size, double con1, double con2)
{
	int j,k;

	for(j=0;j<size;j++)
	{
		for(k=(j+1);k<size;k++)
		{
			out[j*size + k].r += con1 * pcomm[j*size + k].r - con2 * pmat[j*size + k].i;
			out[j*size + k].i += con1 * pcomm[j*size + k].i + con2 * pmat[j*size + k].r;

			out[k*size + j].r = -out[j*size + k].r;
			out[k*size + j].i = out[j*size + k].i;
		}
		out[j*size + j].i += con1 * pcomm[j*size + j].i + con2 * pmat[j*size + j].r;
	}

	return 0;
}

/*
 * adds a hermitian matrix with complex coefficient to output
 * used in 2MM, HMCYM, FreeGaussian, HMC
 */

int AddMat4(doublecomplex *out, doublecomplex *pcomm, int size, double con)
{
	int j,k;

	for(j=0;j<size;j++)
	{
		for(k=(j+1);k<size;k++)
		{
			out[j*size + k].r += -con * pcomm[j*size + k].i;
			out[j*size + k].i += con * pcomm[j*size + k].r;

			out[k*size + j].r = -out[j*size + k].r;
			out[k*size + j].i = out[j*size + k].i;
		}
		out[j*size + j].r += -con * pcomm[j*size + j].i;
	}

	return 0;
}
/*
 * used in Dynamics(8MM), HMC(8MM)
 */
int AddMat4Dyn(doublecomplex *out, doublecomplex *pcomm, int size, double con)
{
	int j,k;

	for(j=0;j<size;j++)
	{
		for(k=(j+1);k<size;k++)
		{
			out[j*size + k].r += -con * pcomm[k*size + j].i;
			out[j*size + k].i += con * pcomm[k*size + j].r;

			out[k*size + j].r = -out[j*size + k].r;
			out[k*size + j].i = out[j*size + k].i;
		}
		out[j*size + j].r += -con * pcomm[j*size + j].i;
	}

	return 0;
}
/*
 * adds hermitian matrix to output
 * used in Dynamics(3MM), Dynamics(8MM), HMCYM, HMC(8MM), HMC(3MM)
 */
int AddMat6(doublecomplex *out, doublecomplex *mat1, double con, int size)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		for(j=(i+1);j<size;j++)
		{
			out[i*size+j].r += con * mat1[i*size+j].r;
			out[i*size+j].i	+= con * mat1[i*size+j].i;

			out[j*size+i].r = out[i*size+j].r;
			out[j*size+i].i = -out[i*size+j].i;
		}
		out[i*size+i].r += con * mat1[i*size+i].r;
	}

	return 0;
}

/*
 * adds hermitian matrix to output
 * BUT adds element (i,j) to element (j,i)!!!
 * used in Dynamics(8MM), HMC(8MM)
 */
int AddMat7(doublecomplex *out, doublecomplex *mat1, double con, int size)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		for(j=(i+1);j<size;j++)
		{
			out[i*size+j].r += con * mat1[j*size+i].r;
			out[i*size+j].i	+= con * mat1[j*size+i].i;

			out[j*size+i].r = out[i*size+j].r;
			out[j*size+i].i = -out[i*size+j].i;
		}
		out[i*size+i].r += con * mat1[i*size+i].r;
	}

	return 0;
}

/*
 * adds matrix hermitian/non-hermitian to output
 */
int AddMat8(doublecomplex *out, doublecomplex *mat1, double con, int size)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			out[i*size+j].r += con * mat1[i*size+j].r;
			out[i*size+j].i	+= con * mat1[i*size+j].i;
		}
	}

	return 0;
}

/*
 * add 2 matrices where the second matrix has a coefficient 'i';
 * used in 2MM, FreeGaussian
 */
int AddMat9(doublecomplex *out, doublecomplex *pmat, doublecomplex *pmat2, int size, double con)
{
	int j,k;

	for(j=0;j<size;j++)
	{
		for(k=0;k<size;k++)
		{
			out[j*size + k].r += con * (pmat[j*size + k].r - pmat2[j*size+k].i);
			out[j*size + k].i += con * (pmat[j*size + k].i + pmat2[j*size+k].r);
		}
	}

	return 0;
}












