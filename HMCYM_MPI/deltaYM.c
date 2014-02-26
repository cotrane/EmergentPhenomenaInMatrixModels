#include<mpi.h>
#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"
#include "MatrixMan.h"

extern MPI_Datatype CMPI_COMPLEX_TYPE;
extern MPI_Op CMPI_ADDCMAT_OP;

int deltaYM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double eps, double mass, int numprocs, int rank)
{
	int i,j,k,a,b,pos1,pos2,numcomm=0;
	double con, tr=0;
	doublecomplex **pcomm, *tmpmat[nummat], *tmpmat2[nummat], **pcommtmp;

	for(i=0;i<nummat;i++)
	{
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmpmat[i], 0, sizeof(doublecomplex)*size*size);
		tmpmat2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmpmat2[i], 0, sizeof(doublecomplex)*size*size);
	}

	for(i=0;i<nummat;i++)
		numcomm += (nummat-1-i);

	pcomm = (doublecomplex**) malloc(sizeof(doublecomplex*)*numcomm);
	pcommtmp = (doublecomplex**) malloc(sizeof(doublecomplex*)*numcomm);
	for(i=0;i<numcomm;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
		pcommtmp[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcommtmp[i], 0, sizeof(doublecomplex)*size*size);
	}

	for(i=0;i<nummat;i++)
	{
		MPI_Bcast(pmat[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);
	}

	for(a=rank; a<(nummat-1); a+=numprocs)							//computes commutator
	{
		pos1=0;
		for(i=0; i<a; i++)
		{
			pos1+=((nummat-1)-i);
		}
		pos2=0;
		for(b=(a+1); b<nummat; b++)
		{
			Comm(pcommtmp[pos1+pos2], pmat[a], pmat[b], size, 1.0);
			pos2++;
		}
	}

	for(i=0;i<numcomm;i++)
	{
		MPI_Reduce(pcommtmp[i], pcomm[i], size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
	}
	for(i=0;i<numcomm;i++)
	{
		MPI_Bcast(pcomm[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);
	}

	for(a=rank;a<nummat;a+=numprocs)								// computes YM Term
	{
		for(b=0;b<nummat;b++)
		{
			if(b!=a)
			{
				if(a<b)
				{
					con=1.0;
					pos1=0;
					for(i=0; i<a; i++)
					{
						pos1+=((nummat-1)-i);
					}
					Commend2(tmpmat2[a], pmat[b], pcomm[pos1+b-1-a], size, eps*size*con);
				}
				if(b<a)
				{
					con=-1.0;
					pos1=0;
					for(i=0; i<b; i++)
					{
						pos1+=((nummat-1)-i);
					}
					Commend2(tmpmat2[a], pmat[b], pcomm[pos1+a-1-b], size, eps*size*con);
				}

			}
		}
	}
	for(i=0;i<nummat;i++)
	{
		MPI_Reduce(tmpmat2[i], tmpmat[i], size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
	}
	if(rank==0)
	{
		for(j=0;j<nummat;j++)
		{
			tr=0;
			for(i=0;i<(size-1);i++)
			{
				(tmpmat[j]+i*size+i)->r -= (tmpmat[j]+(size-1)*size+(size-1))->r;
				tr += (tmpmat[j]+i*size+i)->r;
			}
			(tmpmat[j]+(size-1)*size+(size-1))->r = -tr;
		}
		for(i=0;i<nummat;i++)
		{
			AddMat6(pmom[i], tmpmat[i], 1.0, size);
		}
	}

	for(i=0;i<nummat;i++)
	{
		free(tmpmat[i]);
		free(tmpmat2[i]);
	}
	for(i=0;i<numcomm;i++)
	{
		free(pcomm[i]);
		free(pcommtmp[i]);
	}
	free(pcomm);
	free(pcommtmp);

	return 0;
}

