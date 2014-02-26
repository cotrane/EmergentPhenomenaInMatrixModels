#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include<mpi.h>
#include "MatrixMan.h"


extern MPI_Datatype CMPI_COMPLEX_TYPE;
extern MPI_Op CMPI_ADDCMAT_OP;

int addmat(doublecomplex *pmat[], doublecomplex *pmom[], double eps, int nummat, int size, int numprocs, int rank)
{
	int i,j,k;
	doublecomplex *ptmp[nummat], *ptmp2[nummat];
	double trace=0;

	for(i=0;i<nummat;i++)
	{
		MPI_Bcast(pmom[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);
		ptmp[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(ptmp[i], 0, sizeof(doublecomplex)*size*size);
		ptmp2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(ptmp2[i], 0, sizeof(doublecomplex)*size*size);
	}

	for(i=rank;i<nummat;i+=numprocs)
	{
		trace=0;
		for(j=0;j<(size-1);j++)
		{
			for(k=(j+1);k<size;k++)
			{
				(ptmp[i]+(j*size + k))->r += eps * ((pmom[i]+(k*size + j))->r);
				(ptmp[i]+(j*size + k))->i += eps * (pmom[i]+(k*size + j))->i;

				(ptmp[i]+(k*size + j))->r = (ptmp[i]+(j*size + k))->r;
				(ptmp[i]+(k*size + j))->i = -(ptmp[i]+(j*size + k))->i;
			}
			(ptmp[i]+(j*size + j))->r += eps * (pmom[i]+(j*size + j))->r;
			trace += (ptmp[i]+(j*size + j))->r;
		}
		(ptmp[i]+(size-1)*size + (size-1))->r = -trace;
	}

	for(i=0;i<nummat;i++)
	{
		MPI_Reduce(ptmp[i], ptmp2[i], size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
	}
	if(rank==0)
	{
		for(i=0;i<nummat;i++)
		{
			AddMat6(pmat[i], ptmp2[i], 1.0, size);
		}
	}
	for(i=0;i<nummat;i++)
	{
		free(ptmp[i]);
		free(ptmp2[i]);
	}

	return 0;
}

int addmom(doublecomplex *pmom[], doublecomplex *pmat[], double eps, int nummat, int size, double kappa)
{
	int j,k,l;

	for(j=0;j<nummat;j++)
	{
		for(k=0;k<size;k++)
		{
			for(l=(k+1);l<size;l++)
			{
				(pmom[j]+(k*size + l))->r += -eps * (pmat[j]+(l*size + k))->r;
				(pmom[j]+(k*size + l))->i += -eps * (pmat[j]+(l*size + k))->i;

				(pmom[j]+(l*size + k))->r = (pmom[j]+(k*size + l))->r;
				(pmom[j]+(l*size + k))->i = -(pmom[j]+(k*size + l))->i;
			}
			(pmom[j]+(k*size + k))->r += -eps * (pmat[j]+(k*size + k))->r;
		}
	}
	return 0;
}

double mom(doublecomplex *pmom[], int nummat, int size, double lambda)
{
	int i,j,k;
	double mom=0;

	for(i=0;i<nummat;i++)
	{
		for(j=0;j<(size-1);j++)
		{
			for(k=0;k<size;k++)
			{
				mom += 0.5*((pmom[i]+j*size+k)->r * (pmom[i]+k*size+j)->r - (pmom[i]+j*size+k)->i * (pmom[i]+k*size+j)->i);
			}
		}
		for(k=0;k<(size-1);k++)
		{
			mom += 0.5*((pmom[i]+(size-1)*size+k)->r * (pmom[i]+k*size+(size-1))->r - (pmom[i]+(size-1)*size+k)->i * (pmom[i]+k*size+(size-1))->i);
		}
	}

	return mom;
}
