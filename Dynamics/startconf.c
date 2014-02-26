#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include "RepsSUn.h"
#include "model.h"
#include<clapack.h>
#include<string.h>
#include<math.h>

int gen_randphicplx(doublecomplex *pmat[], int nummat, int size);
int gen_randphicplxtrace(doublecomplex *pmat[], int nummat, int size);

int startconf(doublecomplex *pmat[], double alpha, int nummat, int size, int liegroup, int conf, int *L)
{
	int i,j, N=size*size, Lsu2;
	double a, C2=0;
	doublecomplex *tmpsu3[nummat], *tmpsu2[nummat];

	for(i=0;i<nummat;i++)
	{
		tmpsu2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		tmpsu3[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);

		memset(tmpsu2[i], 0, sizeof(doublecomplex)*N);
		memset(tmpsu3[i], 0, sizeof(doublecomplex)*N);
	}

	/*
	 * for 1MM (Harmonic Oscillator)
	 */
	if(!(liegroup-1))
	{
		gen_randphicplxtrace(pmat, nummat, size);
		return 0;
	}

	a = ((liegroup-2) ? alpha/pow(size, 0.25):alpha);

	*L = RepsSUn(tmpsu3, liegroup, size, nummat, a);
	C2=(liegroup-1.0)/(2.0*liegroup) * *L*(*L+liegroup);
	printf("SU(%d) ground state energy = %f\n", liegroup, -((liegroup-2) ? size:size*size)*liegroup*alpha*alpha*alpha*alpha*C2/12);

	// initialize in SU(3) symmetric configuration; only makes sense for 8MM!!!
	if(conf==0)
	{
		for(i=0;i<nummat;i++)
		{
			memcpy(pmat[i], tmpsu3[i], sizeof(doublecomplex)*N);
		}
	}
	// initialize in random configuration; not very useful for just the Dynamics!
	if(conf==1)
	{
		gen_randphicplx(pmat, nummat, size);
	}
	// initialize in SU(2) symmetric configuration
	if(conf==2)
	{
		for(i=0;i<nummat;i++)
		{
			memset(pmat[i], 0, sizeof(doublecomplex)*N);
		}
		Lsu2=RepsSUn(pmat, 2, size, 3, a);

		C2=1.0/4.0 * Lsu2*(Lsu2+2);
		printf("SU(2) ground state energy = %f\n", -size*size*2*a*a*a*a*C2/12);
	}

	// ---------------------------------- for 8MM only !!! ----------------------------------------------
	// see thesis for different configurations for SU(2) symmetric state in 8MM
	// initialize in SU(2) symmetric configuration formed by pmat[3], pmat[4], and pmat[2] and pmat[7]
	if(conf==3)
	{
		Lsu2=RepsSUn(tmpsu2, 2, size, 3, a);

		C2=1.0/4.0 * Lsu2*(Lsu2+2);
		printf("SU(2) ground state energy = %f\n", -size*size*2*a*a*a*a*C2/12);

		memset(pmat[3], 0, sizeof(doublecomplex)*N);
		memset(pmat[4], 0, sizeof(doublecomplex)*N);

		memcpy(pmat[3], tmpsu2[0], sizeof(doublecomplex)*N);
		memcpy(pmat[4], tmpsu2[1], sizeof(doublecomplex)*N);
		memset(pmat[0], 0, sizeof(doublecomplex)*N);
		memset(pmat[1], 0, sizeof(doublecomplex)*N);
		memset(pmat[5], 0, sizeof(doublecomplex)*N);
		memset(pmat[6], 0, sizeof(doublecomplex)*N);

		for(i=0;i<size;i++)
		{
			for(j=0;j<size;j++)
			{
				(pmat[2]+(i*size+j))->r = 0.5*(tmpsu2[2]+(i*size+j))->r;

				(pmat[7]+(i*size+j))->r = sqrt(3.0)*0.5*(tmpsu2[2]+(i*size+j))->r;
			}
		}
	}

	// initialize in SU(2) symmetric configuration formed by pmat[5], pmat[6], and pmat[2] and pmat[7]
	if(conf==4)
	{
		Lsu2=RepsSUn(tmpsu2, 2, size, 3, a);

		C2=1.0/4.0 * Lsu2*(Lsu2+2);
		printf("SU(2) ground state energy = %f\n", -size*size*2*a*a*a*a*C2/12);

		memset(pmat[5], 0, sizeof(doublecomplex)*N);
		memset(pmat[6], 0, sizeof(doublecomplex)*N);

		memcpy(pmat[5], tmpsu2[0], sizeof(doublecomplex)*N);
		memcpy(pmat[6], tmpsu2[1], sizeof(doublecomplex)*N);
		memset(pmat[0], 0, sizeof(doublecomplex)*N);
		memset(pmat[1], 0, sizeof(doublecomplex)*N);
		memset(pmat[3], 0, sizeof(doublecomplex)*N);
		memset(pmat[4], 0, sizeof(doublecomplex)*N);

		for(i=0;i<size;i++)
		{
			for(j=0;j<size;j++)
			{
				(pmat[2]+(i*size+j))->r = -0.5*(tmpsu2[2]+(i*size+j))->r;

				(pmat[7]+(i*size+j))->r = sqrt(3.0)*0.5*(tmpsu2[2]+(i*size+j))->r;
			}
		}
	}

	return 0;
}
