#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "RepsSUn.h"
#include "model.h"
#include"hmc.h"

int startconf(doublecomplex *pmat[], double alpha, int nummat, int size, int liegroup, int numrep, int *rep, int conf, int *L, int diag,
		double prop1, double prop2)
{
	int i,j,k,l, N=size*size, Lsu2, repsize=size, repsizeold=0;
	double a, C2=0;
	doublecomplex *tmpsu3[nummat], *tmpsu2[nummat];

	for(i=0;i<nummat;i++)
	{
		tmpsu2[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);
		tmpsu3[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*N);

		memset(tmpsu2[i], 0, sizeof(doublecomplex)*N);
		memset(tmpsu3[i], 0, sizeof(doublecomplex)*N);
	}

	// determine value of constant 'a' depending on model
	a = ((liegroup-2) ? alpha/pow(size, 0.25):alpha);
	*L = WeylDimFormula(size, liegroup);

	/*
	 * conf==0:
	 * if numrep==0 -> start with irrep of SU(liegroup) of size N
	 * if numrep!=0 -> start with irrep of SU(liegroup) of size specified by numrep;
	 * 	numrep defined as number of boxes in Yang-diagram for SU(liegroup)
	 * 	rep is put into bigger matrix starting from element (0,0) to (repsize,repsize);
	 */
	if(conf==0)
	{
		if(numrep==1 && rep[0] == 0)
		{
			for(i=0;i<nummat;i++)
			{
				RepsSUn(pmat, liegroup, size, nummat, a);
			}
		}
		else
		{
			for(i=0;i<numrep;i++)
			{
				for(j=0;j<nummat;j++)
				{
					free(tmpsu3[j]);
				}
				repsize = (rep[i]+1)*(rep[i]+2)/2.0;
				for(j=0;j<nummat;j++)
				{
					tmpsu3[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*repsize*repsize);
					memset(tmpsu3[j], 0, sizeof(doublecomplex)*repsize*repsize);
				}
				RepsSUn(tmpsu3, liegroup, repsize, nummat, a);
				for(l=0;l<nummat;l++)
				{
					for(j=0;j<repsize;j++)
					{
						for(k=0;k<repsize;k++)
						{
							(pmat[l]+(repsizeold*size+repsizeold + j*size+k))->r = (tmpsu3[l]+j*repsize+k)->r;
							(pmat[l]+(repsizeold*size+repsizeold + j*size+k))->i = (tmpsu3[l]+j*repsize+k)->i;
						}
					}
				}
				repsizeold+=repsize;
			}
		}
		C2=(liegroup-1.0)/(2.0*liegroup) * *L*(*L+liegroup);
		printf("SU(3) ground state = %f\n",-size*liegroup*alpha*alpha*alpha*alpha*C2/12);
	}
	/*
	 * conf==1: generate random starting configuration
	 */
	if(conf==1)
	{
		gen_randphicplx(pmat, size);
	}
	/*
	 * cold start in SU(2) rep; size of SU(2) embedding is specified by numrep;
	 * if repsize<size: embedding is put into submatrix (0,0) to (repsize,repsize)
	 * for 8MM there are 3 different ways of embedding SU(2); in matrix 1-2-3, 4-5-(3,8), 6-7-(3,8)
	 * conf==2 uses first embedding
	 * conf==3 		second
	 * conf==4		third
	 */
	if(conf==2)
	{
		if(rep[0]==0)
		{
			for(i=0;i<nummat;i++)
			{
				memset(pmat[i], 0, sizeof(doublecomplex)*N);
			}
			Lsu2=RepsSUn(pmat, 2, size, 3, a);
			repsize=size;
		}
		else
		{
			for(i=0;i<numrep;i++)
			{
				for(j=0;j<nummat;j++)
				{
					free(tmpsu2[j]);
				}
				repsize=1;
				for(j=1;j<2;j++)
				{
					repsize = repsize*(rep[i]+j)/(j*1.0);
				}
				for(j=0;j<nummat;j++)
				{
					tmpsu2[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*repsize*repsize);
					memset(tmpsu2[j], 0, sizeof(doublecomplex)*repsize*repsize);
				}
				RepsSUn(tmpsu2, 2, repsize, 3, a);
				for(l=0;l<3;l++)
				{
					for(j=0;j<repsize;j++)
					{
						for(k=0;k<repsize;k++)
						{
							(pmat[l]+(repsizeold*size+repsizeold + j*size+k))->r = prop1*(tmpsu2[l]+j*repsize+k)->r;
							(pmat[l]+(repsizeold*size+repsizeold + j*size+k))->i = prop1*(tmpsu2[l]+j*repsize+k)->i;
						}
					}
				}
				repsizeold+=repsize;
			}
			Lsu2 = WeylDimFormula(repsize, 2);
			printf("Lsu2=%d\n", Lsu2);
//			Lsu2=RepsSUn(tmpsu3, 2, repsize, 3, a);
		}
		C2=1.0/4.0 * Lsu2*(Lsu2+2);
		printf("SU(2) ground state = %f\n", -size*repsize*2*a*a*a*a*C2/12);
	}
	if(conf==3)
	{
		Lsu2=RepsSUn(tmpsu2, 2, size, 3, a);

		C2=1.0/4.0 * Lsu2*(Lsu2+2);
		printf("SU(2) ground state = %f\n", -size*size*2*a*a*a*a*C2/12);

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
	if(conf==4)
	{
		Lsu2=RepsSUn(tmpsu2, 2, size, 3, a);

		C2=1.0/4.0 * Lsu2*(Lsu2+2);
		printf("SU(2) ground state = %f\n", -size*size*2*a*a*a*a*C2/12);

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
	/*
	 * cold start in SU(2)xU(1) rep
	 * it is necessary to use conf==2 as U(1) factor is always put into 8th matrix!!!
	 * size of embedding is specified over repsize-> see conf==2
	 */
	if(diag!=0)
	{
		memset(pmat[7], 0, sizeof(doublecomplex)*size*size);
		for(i=0;i<repsize;i++)
		{
			(pmat[7]+i*size+i)->r = a*prop2*1.0;
		}
		for(i=repsize;i<size;i++)
		{
			(pmat[7]+i*size+i)->r = -a*prop2*repsize/(size*1.0-repsize*1.0)*1.0;
		}
	}

	return 0;
}
