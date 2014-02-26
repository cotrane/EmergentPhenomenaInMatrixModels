#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "prng.h"
#include "montecarlo.h"
#include "hmc.h"
#include<mpi.h>
#include "MatrixMan.h"

void CMPI_addmat (doublecomplex *a, doublecomplex *b, int *len, MPI_Datatype *type);
extern MPI_Datatype CMPI_COMPLEX_TYPE;
extern MPI_Op CMPI_ADDCMAT_OP;

int deltaX8MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps, double C2, double C3, double lambda, int numprocs, int rank)
{
	int i,j,a,b,c,save, savea=0, saveb=0, savec=0, pos1, pos2, turn=0;
	double mat_c[36], con, tr;
	doublecomplex *pcomm[28], *pcommtmp[28], *tmpmat[nummat], *ptmp[nummat], *ptmp2[nummat];

	for(i=0;i<nummat;i++)
	{
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmpmat[i], 0, sizeof(doublecomplex)*size*size);
		ptmp[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(ptmp[i], 0, sizeof(doublecomplex)*size*size);
		ptmp2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(ptmp2[i], 0, sizeof(doublecomplex)*size*size);
	}
	for(i=0;i<28;i++)
	{
		pcommtmp[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcommtmp[i], 0, sizeof(doublecomplex)*size*size);
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	for(i=0;i<nummat;i++)
	{
		MPI_Bcast(pmat[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);
	}

	if(lambda!=0)
	{
		deltaS1(pmat, tmpmat, nummat, size, eps, C2, C3, lambda, coupling, numprocs, rank);
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

	for(i=0;i<28;i++)
	{
		MPI_Reduce(pcommtmp[i], pcomm[i], size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
	}
	for(i=0;i<28;i++)
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
					Commend2(ptmp[a], pmat[b], pcomm[pos1+b-1-a], size, eps*size*con);
				}
				else if(b<a)
				{
					con=-1.0;
					pos1=0;
					for(i=0; i<b; i++)
					{
						pos1+=((nummat-1)-i);
					}
					Commend2(ptmp[a], pmat[b], pcomm[pos1+a-1-b], size, eps*size*con);
				}

			}
		}
	}

	structconst(mat_c);
	for(i=4*rank;i<=32;i+=4*numprocs)
	{
		a=mat_c[i];
		b=mat_c[i+1];
		c=mat_c[i+2];
		con=mat_c[i+3];

		savea = a;
		//abc
		pos1=0;
		for(j=0;j<b;j++)
		{
			pos1+=((nummat-1)-j);
		}
		AddMat4Dyn(ptmp[a], pcomm[pos1+c-1-b], size, -2.0*eps*coupling*con*size/pow(size, 0.25));

		//bca
		save = a;
		a=b;
		b=c;
		c=save;

		if(b>c)
		{
			con = -con;
			saveb=b;
			savec=c;
			b=c;
			c=saveb;
			turn=1;
		}

		pos1=0;
		for(j=0;j<b;j++)
		{
			pos1+=((nummat-1)-j);
		}
		AddMat4Dyn(ptmp[a], pcomm[pos1+c-1-b], size, -2.0*eps*coupling*con*size/pow(size, 0.25));

		if(turn==1)
		{
			b=saveb;
			c=savec;
			con=-con;
			turn=0;
		}

		//cab
		save = a;
		a=b;
		b=c;
		c=save;

		if(b>c)
		{
			con = -con;
			saveb=b;
			savec=c;
			b=c;
			c=saveb;
			turn=1;
		}

		pos1=0;
		for(j=0;j<b;j++)
		{
			pos1+=((nummat-1)-j);
		}
		AddMat4Dyn(ptmp[a], pcomm[pos1+c-1-b], size, -2.0*eps*coupling*con*size/pow(size, 0.25));
	}

	for(i=0;i<nummat;i++)
	{
		MPI_Reduce(ptmp[i], ptmp2[i], size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
	}

	if(rank==0)
	{
		for(i=0;i<nummat;i++)
		{
			AddMat6(tmpmat[i], ptmp2[i], 1.0, size);
		}
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
		free(ptmp[i]);
		free(ptmp2[i]);
	}
	for(i=0;i<28;i++)
	{
		free(pcomm[i]);
		free(pcommtmp[i]);
	}

	return 0;
}
