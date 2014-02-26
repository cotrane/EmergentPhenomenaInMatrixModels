/*
 * genD.c
 *
 *  Created on: 8 Feb 2013
 *      Author: tkaltenbrunner
 */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<f2c.h>
#include<clapack.h>
#include"MatrixMan.h"
#include"montecarlo.h"
#include<mpi.h>

/*
 * computes Dirac operator D = \gamma_{\mu} \otimes [ X_{\mu}, . ]; D = 4*N^2*N^2 matrix
 */

extern MPI_Datatype CMPI_COMPLEX_TYPE;
extern MPI_Op CMPI_ADDCMAT_OP;

int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat, int numprocs, int rank)
{
	int i,j,k,l,m,n,I,J,mu;
	int sizeD = size*size;
	doublecomplex *ptmp, *ptmp2;

	for(i=0;i<nummat;i++)
	{
		MPI_Bcast(pmat[i], size*size, CMPI_COMPLEX_TYPE, 0, MPI_COMM_WORLD);
	}
	ptmp = (doublecomplex*) malloc(sizeof(doublecomplex)*sizegamma*sizegamma*sizeD*sizeD);
	memset(ptmp, 0, sizeof(doublecomplex)*sizegamma*sizegamma*sizeD*sizeD);
	ptmp2 = (doublecomplex*) malloc(sizeof(doublecomplex)*sizegamma*sizegamma*sizeD*sizeD);
	memset(ptmp2, 0, sizeof(doublecomplex)*sizegamma*sizegamma*sizeD*sizeD);

	//dirac op as 2*N2x2*N2 matrix
	for(mu=rank;mu<nummat;mu+=numprocs){
		for(i=0;i<size;i++){
			for(j=0;j<size;j++){
				for(k=0;k<size;k++){
					if(i==j && j==k)
						continue;
					for(m=0;m<sizegamma;m++){
						for(n=0;n<sizegamma;n++){
							I=i*size+j;J=k*size+j;
							ptmp[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].r += (gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+k*size+i)->r - (gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+k*size+i)->i;
							ptmp[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].i += (gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+k*size+i)->r + (gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+k*size+i)->i;

							I=i*size+j;J=i*size+k;
							ptmp[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].r += -((gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+j*size+k)->r - (gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+j*size+k)->i);
							ptmp[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].i += -((gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+j*size+k)->r + (gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+j*size+k)->i);
						}
					}
				}
			}
		}
	}
	MPI_Reduce(ptmp, ptmp2, sizegamma*sizegamma*sizeD*sizeD, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);

	if(rank==0)
		memcpy(out, ptmp2, sizeof(doublecomplex)*sizegamma*sizegamma*sizeD*sizeD);

	free(ptmp);
	free(ptmp2);

//------------------		testing		---------------------------------------------
//	int l;
//	double trr=0, tri=0;
//	doublecomplex *tmp2, *tmp3;
//	tmp2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
//	memset(tmp2, 0, sizeof(doublecomplex)*size*size);
//	tmp3 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
//	memset(tmp3, 0, sizeof(doublecomplex)*size*size);
//
//	memcpy(tmp3, pmat[1], sizeof(doublecomplex)*size*size);
////	AddMat5(tmp3, pmat[1], pmat[2],1,1,size);
//	printmat(tmp3, size);
//
//	for(i=0;i<size;i++){
//		for(j=0;j<size;j++){
//			for(k=0;k<size;k++){
//				for(l=0;l<size;l++){
//					I=i*size+j;J=k*size+l;
//					if((i==(size-1) && j==(size-1)) || (k==(size-1) && l==(size-1)))
//						continue;
//					tmp2[i*size+j].r += out[I*sizeD + J].r * tmp3[k*size+l].r - out[I*sizeD + J].i * tmp3[k*size+l].i;
//					tmp2[i*size+j].i += out[I*sizeD + J].r * tmp3[k*size+l].i + out[I*sizeD + J].i * tmp3[k*size+l].r;
//					if(i==j){
//						trr += tmp2[i*size+j].r;
//						tri += tmp2[i*size+j].i;
//					}
//				}
//			}
//		}
//	}
//	tmp2[(size-1)*size+(size-1)].r = -trr;
//	tmp2[(size-1)*size+(size-1)].i = -tri;
//
//	printmat(tmp2, size);

//	for(mu=0;mu<nummat;mu++)
//		free(tmp[mu]);
//	free(tmp);

	return 0;
}
