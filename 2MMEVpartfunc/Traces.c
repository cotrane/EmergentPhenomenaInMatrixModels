/*
 * Traces.c
 *
 *  Created on: 9 Jan 2013
 *      Author: tkaltenbrunner
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "MatrixMan.h"

extern FILE *outTraces;

/*
 * compute trace observables
 */
int Traces(doublecomplex **pmat, int size, int nummat, int loopnumber, char *pathout)
{
	int i;
	static int init=0;
	static doublecomplex *pmatX2, *pmatY2, *pmatX4, *pmatY4, *pmatX6, *pmatY6, *pmatX8, *pmatY8;
	static doublecomplex *pmatXY, *pmatX2Y2, *pmatX4Y2, *pmatX2Y4, *pmatX4Y4, *pmatX6Y2, *pmatX6Y4, *pmatX6Y6, *pmatX8Y2, *pmatX8Y4, *pmatX8Y6, *pmatX8Y8;
	double trX2=0, trXi2=0, trXY=0, trX2Y2=0, trXYXY=0, trX2Y=0, trX2Y4=0, trX4Y2=0, trX4Y4=0;

	if(init==0)
	{
		pmatX2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatY2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatY4 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX4 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatY6 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX6 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatY8 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX8 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatXY = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX2Y2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX4Y2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX2Y4 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX4Y4 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX6Y2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX6Y4 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX6Y6 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX8Y2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX8Y4 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX8Y6 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		pmatX8Y8 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);

		init=1;
	}

	memset(pmatX2, 0, sizeof(doublecomplex)*size*size);
	memset(pmatY2, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX4, 0, sizeof(doublecomplex)*size*size);
	memset(pmatY4, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX6, 0, sizeof(doublecomplex)*size*size);
	memset(pmatY6, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX8, 0, sizeof(doublecomplex)*size*size);
	memset(pmatY8, 0, sizeof(doublecomplex)*size*size);
	memset(pmatXY, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX2Y2, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX4Y2, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX2Y4, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX4Y4, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX6Y2, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX6Y4, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX6Y6, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX8Y2, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX8Y4, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX8Y6, 0, sizeof(doublecomplex)*size*size);
	memset(pmatX8Y8, 0, sizeof(doublecomplex)*size*size);

	Multi5(pmatX2, pmat[0], pmat[0], size, 1.0);
	Multi5(pmatY2, pmat[1], pmat[1], size, 1.0);
	Multi5(pmatX4, pmatX2, pmatX2, size, 1.0);
	Multi5(pmatY4, pmatY2, pmatY2, size, 1.0);
	Multi5(pmatX6, pmatX4, pmatX2, size, 1.0);
	Multi5(pmatY6, pmatY4, pmatY2, size, 1.0);
	Multi5(pmatX8, pmatX4, pmatX4, size, 1.0);
	Multi5(pmatY8, pmatY4, pmatY4, size, 1.0);
	Multi5(pmatXY, pmat[0], pmat[1], size, 1.0);
	Multi5(pmatX2Y2, pmatX2, pmatY2, size, 1.0);
	Multi5(pmatX4Y2, pmatX4, pmatY2, size, 1.0);
	Multi5(pmatX2Y4, pmatX2, pmatY4, size, 1.0);
	Multi5(pmatX4Y4, pmatX4, pmatY4, size, 1.0);
	Multi5(pmatX6Y2, pmatX6, pmatY2, size, 1.0);
	Multi5(pmatX6Y4, pmatX6, pmatY4, size, 1.0);
	Multi5(pmatX6Y6, pmatX6, pmatY6, size, 1.0);
	Multi5(pmatX8Y2, pmatX8, pmatY2, size, 1.0);
	Multi5(pmatX8Y4, pmatX8, pmatY4, size, 1.0);
	Multi5(pmatX8Y6, pmatX8, pmatY6, size, 1.0);
	Multi5(pmatX8Y8, pmatX8, pmatY8, size, 1.0);

	for(i=0;i<nummat;i++)
	{
		diagMulti(&trXi2, pmat[i], pmat[i], size, 1.0);
	}
	diagMulti(&trX2, pmat[0], pmat[0], size, 1.0);
	diagMulti(&trXY, pmat[0], pmat[1], size, 1.0);
	diagMulti(&trXYXY, pmatXY, pmatXY, size, 1.0);
	diagMulti(&trX2Y2, pmatX2, pmatY2, size, 1.0);
	diagMulti(&trX2Y, pmatX2, pmat[1], size, 1.0);
	diagMulti(&trX4Y2, pmatX4, pmatY2, size, 1.0);
	diagMulti(&trX2Y4, pmatX2, pmatY4, size, 1.0);
	diagMulti(&trX4Y4, pmatX4, pmatY4, size, 1.0);

	fprintf(outTraces, "%f\t%f\t%f\t%f\t%f\t%f\n",  trX2/size, trXi2/size, trXY/size, trXYXY/size, trX2Y2/size, trX2Y/size);
	fprintf(outTraces, "%f\t%f\t%f\t%f\n",  trX2Y2/size, trX2Y4/size, trX4Y2/size, trX4Y4/size);

	fflush(outTraces);

	return 0;
}
