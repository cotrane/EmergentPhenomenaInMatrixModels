#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>
#include"MatrixMan.h"

int Ops(int **ppopladder, int N);
int GenStates(int **ppstate, int **ppopladder, int **ppMultiplet, int N, int dim);
int GenMultiplets(int **ppMultiplet, int **ppMultinum, int N, int dim, int numMultis);
int GenDiagonal(doublecomplex **ppmat, double **pplambda, int **ppstate, int **ppMultiplet, int **ppMultinum, int N, int dim, double a);
int GenOffDiagonal(doublecomplex **ppmat, int **ppMultiplet, int **ppMultinum, int N, int dim, double a);
int WeylDimFormula(int size, int liegroup);

/*
 * generates totally symmetric representations of SU(N):
 * ppmat ... output matrices
 * N ... liegroup
 * Matrix_Size ... size of representation
 * Nummat ... number of matrices; must correspond to liegroup!
 * a ... coefficient of matrices; 1.0 for representation; 1/2.0 for generators
 */
int RepsSUn(doublecomplex *ppmat[], int N, int Matrix_Size, int Nummat, double a)
{
	int i;
	int **ppstate, *rep, **ppMultiplet, **ppopladder, **ppMultinum;
	int dim=0, numMultis=0;
	double **pplambda;
	double C2=0, C2analyt=0;
	doublecomplex *pcomm[3], trcomm[3];


	rep = (int*) malloc(sizeof(int)*(N-1));
	memset(rep, 0, sizeof(int)*(N-1));

	rep[0] = WeylDimFormula(Matrix_Size, N);
	dim = Matrix_Size;

	numMultis=100;

	pplambda = (double**) malloc(sizeof(double*)*(N-1));
	for(i=0;i<(N-1);i++)
	{
		pplambda[i] = (double*) malloc(sizeof(double)*N);
	}

	ppMultiplet = (int**) malloc(sizeof(int*)*((N*N-1)-(N-1)));
	ppMultinum = (int**) malloc(sizeof(int*)*((N*N-1)-(N-1)));
	ppopladder = (int**) malloc(sizeof(int*)*((N*N-1)-(N-1)));
	for(i=0;i<((N*N-1)-(N-1));i++)
	{
		ppMultiplet[i] = (int*) malloc(sizeof(int)*dim);
		ppMultinum[i] = (int*) malloc(sizeof(int)*numMultis);
		ppopladder[i] = (int*) malloc(sizeof(int)*(N));
	}

	ppstate = (int**) malloc(sizeof(int*)*dim); // not ppstate[0][1];, but pTmp = ppstate[0]' pTmp[1] = ...
	for(i=0;i<dim;i++)
	{
		ppstate[i] = (int*) malloc(sizeof(int)*(N));
	}

	for(i=0;i<3;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*dim*dim);
	}


	for(i=0;i<(((N*N-1)-(N-1)));i++)
	{
		memset(ppMultinum[i], 0, sizeof(int)*numMultis);
		memset(ppMultiplet[i], 0, sizeof(int)*dim);
		memset(ppopladder[i], 0, sizeof(int)*(N));
	}
	for(i=0;i<(N-1);i++)
	{
		memset(pplambda[i], 0, sizeof(double)*N);
	}
	for(i=0;i<dim;i++)
	{
		memset(ppstate[i], 0, sizeof(int)*(N));
	}
	for(i=0;i<3;i++)
	{
		memset(pcomm[i], 0, dim*dim*sizeof(doublecomplex));
	}
	memset(trcomm,0,3*sizeof(doublecomplex));


	Ops(ppopladder, N);				//generates operators

	for(i=0;i<(N-1);i++)			//highest weight state
	{
		*(ppstate[0]+i) = rep[i];
	}

	GenStates(ppstate, ppopladder, ppMultiplet, N, dim);

	GenMultiplets(ppMultiplet, ppMultinum, N, dim, numMultis);

	GenOffDiagonal(ppmat, ppMultiplet, ppMultinum, N, dim, a);

	GenDiagonal(ppmat, pplambda, ppstate, ppMultiplet, ppMultinum, N, dim, a);

	// ----------------- for testing purposes -------------------------------------------

	for(i=0;i<Nummat;i++)
	{
		diagMulti(&C2, ppmat[i], ppmat[i], dim, 1.0);
	}
	C2analyt = a*a*dim/2.0 *(1.0*N*rep[0]+rep[0]*(rep[0]+1.0-2.0) - 1.0/N *rep[0]*rep[0]);


	if((C2-C2analyt)>0.00000001)
	{
		printf("\n----------------------------------------------------------------------------\n");
		printf("test for RepsSUn: X2:= (X_a)^2 = a^2 C^{adj} C_2 Tr{1} / 2.0\n");
		printf("X2=%lf\tX2analyt=%lf\n", C2, C2analyt);
		printf("error in the computation of the matrices! C2 != C2analyt\n");
		printf("----------------------------------------------------------------------------\n");
	}

	return rep[0];
}
