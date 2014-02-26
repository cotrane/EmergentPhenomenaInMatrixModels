/*
 * Inverse.c
 *
 *  Created on: 31 Jan 2013
 *      Author: tkaltenbrunner
 */


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<f2c.h>
#include<clapack.h>

int printmat(char *use, doublecomplex *mat, int size);

int Inverse_posDef(doublecomplex *out, doublecomplex *in, integer size)
{
	int i,j;
	char uplo='U';
	integer info=0;

	memcpy(out, in, sizeof(doublecomplex)*size*size);

	zpotrf_(&uplo, &size, out, &size, &info);
	if(info!=0){
		printf("Error in computing Cholesky factorization! Info=%ld\n", (long int) info);
		info=0;
	}
	zpotri_(&uplo, &size, out, &size, &info);
	if(info!=0){
		printf("Error in computing the Inverse! Info=%ld\n", (long int) info);
	}

	for(i=0;i<size;i++){
		for(j=(i+1);j<size;j++){
			out[i*size+j].r = out[j*size+i].r;
			out[i*size+j].i = -out[j*size+i].i;
		}
	}

	return 0;
}

int Inverse_indef(doublecomplex *out, doublecomplex *in, integer size)
{
	int i,j;
	char uplo='U';
	integer info=0, *ipiv;
	doublecomplex *work;
	integer lwork = 64*size; //???

	ipiv = (integer*) malloc(sizeof(integer)*size);
	work = (doublecomplex*) malloc(sizeof(doublecomplex)*lwork);
	memcpy(out, in, sizeof(doublecomplex)*size*size);

	zhetrf_(&uplo, &size, out, &size, ipiv, work, &lwork, &info);
	if(info!=0){
		printf("Error in computing Cholesky factorization! Info=%ld\n", (long int) info);
		info=0;
	}
	zhetri_(&uplo, &size, out, &size, ipiv, work, &info);
	if(info!=0){
		printf("Error in computing the Inverse! Info=%ld\n", (long int) info);
	}

	for(i=0;i<size;i++){
		for(j=(i+1);j<size;j++){
			out[i*size+j].r = out[j*size+i].r;
			out[i*size+j].i = -out[j*size+i].i;
		}
	}

	return 0;
}
