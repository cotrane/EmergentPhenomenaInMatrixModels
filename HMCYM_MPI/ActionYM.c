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

double actionYM(doublecomplex *pmat[], int size, int nummat, double mass, int numprocs, int rank)
{
	int i, numcomm=0,pos1,pos2,mat1,mat2;
	double trYM=0, tmpYM=0;
	doublecomplex **pcomm;

	// broadcast matrices pmat to all different nodes
	for(i=0;i<nummat;i++)
	{
		MPI_Bcast(pmat[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);
	}

	for(i=0;i<nummat;i++)
		numcomm += (nummat-1-i);

	pcomm = (doublecomplex**) malloc(sizeof(doublecomplex*)*numcomm);
	for(i=0;i<numcomm;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	for(mat1=rank; mat1<(nummat-1); mat1+=numprocs)
	{
		pos1=0;
		for(i=0; i<mat1; i++)
		{
			pos1+=((nummat-1)-i);
		}
		pos2=0;
		for(mat2=(mat1+1); mat2<nummat; mat2++)
		{
			Comm(pcomm[pos1+pos2], pmat[mat1], pmat[mat2], size, 1.0);
			diagMulti(&tmpYM, pcomm[pos1+pos2], pcomm[pos1+pos2], size, -size/2.0);
			pos2++;
		}
	}
	// collect results from computation by each node in master-node rank = 0
	MPI_Reduce(&tmpYM, &trYM, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	for(i=0;i<numcomm;i++)
	{
		free(pcomm[i]);
	}
	free(pcomm);

	return trYM;
}


