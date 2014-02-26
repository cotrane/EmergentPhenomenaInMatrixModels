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

/*
 * computes Dirac operator D = \gamma_{\mu} \otimes [ X_{\mu}, . ]; D = 4*N^2*N^2 matrix
 */
int genD(doublecomplex *out, doublecomplex **pmat, doublecomplex **gamma, int size, int sizegamma, int nummat)
{
	int i,j,k,m,n,I,J,mu;
	int sizeD = size*size;

	//dirac op as N2xN2 matrix
	for(mu=0;mu<nummat;mu++){
		for(i=0;i<size;i++){
			for(j=0;j<size;j++){
				for(k=0;k<size;k++){
					if(i==j && j==k)
						continue;
					for(m=0;m<sizegamma;m++){
						for(n=0;n<sizegamma;n++){
							I=i*size+j;J=k*size+j;
							out[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].r += (gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+k*size+i)->r - (gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+k*size+i)->i;
							out[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].i += (gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+k*size+i)->r + (gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+k*size+i)->i;

							I=i*size+j;J=i*size+k;
							out[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].r += -((gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+j*size+k)->r - (gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+j*size+k)->i);
							out[(I*sizegamma+m)*(sizegamma*sizeD) + (J*sizegamma+n)].i += -((gamma[mu]+m*sizegamma+n)->i * (pmat[mu]+j*size+k)->r + (gamma[mu]+m*sizegamma+n)->r * (pmat[mu]+j*size+k)->i);
						}
					}
				}
			}
		}
	}

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
