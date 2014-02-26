/*
 * ConstructAdjointAction2.c
 *
 *  Created on: 2 Nov 2012
 *      Author: tkaltenbrunner
 */
/*
 * ConstructAdjointRep.c
 *
 *  Created on: 29 Oct 2012
 *      Author: tkaltenbrunner
 */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<f2c.h>
#include "MatrixMan.h"
#include "montecarlo.h"
#include "RepsSUn.h"
#include "hmc.h"


doublecomplex **basis;
structconst *StructConst;

int ConstructAdjointAction(doublecomplex **out, doublecomplex **pmat, int size, int nummat)
{
	int i,j,k, numbasis = size*size-1;
	int a,b,c;
	static int init=0;
	double trace;
	doublecomplex *comm, *test;

	comm = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
	test = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);

	if(init==0)
	{
		StructConst = (structconst*) malloc(sizeof(structconst)*3);
		basis = (doublecomplex**) malloc(sizeof(doublecomplex*)*(numbasis));
		for(i=0;i<(numbasis);i++)
		{
			basis[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
			memset(basis[i],0, sizeof(doublecomplex)*size*size);
		}
		StructConst[0].numEntries = 3;
		RepsSUn(basis, size, size, numbasis, 2.0);
//		for(i=0;i<numbasis;i++)
//			printmat("basis", basis[i], size);
		findStructConst(&StructConst, size);
//		printf("numEntries=%d\n", (StructConst)[0].numEntries);
//		for(i=0;i<((StructConst)[0].numEntries);i++)
//			printf("%d-%d-%d: %f\n", (StructConst)[i].a, (StructConst)[i].b, (StructConst)[i].c, (StructConst)[i].val);
		init=1;
	}

	for(i=0;i<(StructConst[0].numEntries);i++)
	{
		a=StructConst[i].a;
		b=StructConst[i].b;
		c=StructConst[i].c;
		for(j=0;j<nummat;j++)
		{
			trace=0;
			diagMultiI(&trace, pmat[j], basis[a], size, StructConst[i].val);
			(out[j]+c*numbasis+b)->i += -trace;
			(out[j]+b*numbasis+c)->i += trace;

			trace=0;
			diagMultiI(&trace, pmat[j], basis[b], size, StructConst[i].val);
			(out[j]+a*numbasis+c)->i += -trace;
			(out[j]+c*numbasis+a)->i += trace;

			trace=0;
			diagMultiI(&trace, pmat[j], basis[c], size, StructConst[i].val);
			(out[j]+b*numbasis+a)->i += -trace;
			(out[j]+a*numbasis+b)->i += trace;
		}
	}

	free(comm);
	free(test);

	return 0;
}

