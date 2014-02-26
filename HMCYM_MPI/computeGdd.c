/*
 * computeGdd.c
 *
 *  Created on: 22 Oct 2013
 *      Author: tkaltenbrunner
 */

#include<mpi.h>
#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include"MatrixMan.h"
#include<string.h>

MPI_Datatype	CMPI_COMPLEX_TYPE;
MPI_Op			CMPI_ADDCMAT_OP;

int computeGdd(doublecomplex *out, doublecomplex **in, int nummat, int size, int rank, int numprocs)
{
	int i,j;
	doublecomplex val, *tmp;

	tmp = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
	memset(tmp, 0, sizeof(doublecomplex)*size*size);

	for(i=0;i<nummat;i++)
		MPI_Bcast(in[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);

	for(i=rank;i<nummat;i+=numprocs){
		for(j=i;j<nummat;j++){
			val.r=0; val.i=0;
			diagMultiCplx(&val, in[i], in[j], size, 1.0);
			if(i!=j){
				tmp[i*nummat+j].r = val.r;
				tmp[i*nummat+j].i = val.i;

				tmp[j*nummat+i].r = val.r;
				tmp[j*nummat+i].i = -val.i;
			}
			else{
				tmp[i*nummat+i].r = val.r;
				tmp[i*nummat+i].i = val.i;
			}
		}
	}

	MPI_Reduce(tmp, out, size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);

	return 0;
}
