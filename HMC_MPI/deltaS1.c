#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include<mpi.h>
#include "MatrixMan.h"

void CMPI_addmat (doublecomplex *a, doublecomplex *b, int *len, MPI_Datatype *type);
extern MPI_Datatype CMPI_COMPLEX_TYPE;
extern MPI_Op CMPI_ADDCMAT_OP;

int deltaS1(doublecomplex *pmat[], doublecomplex *tmpmat[], int nummat, int size, double eps, double C2, double C3, double lambda, double alpha, int numprocs, int rank)
{
	int i,j,k,a,b,c,e,f,g;
	double mat_c2[32], mat_c22[64];
	double con=0, con2=0;
	doublecomplex *tmpmat1, *tmpmat2, *tmp[nummat];

	tmpmat1 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
	tmpmat2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
	for(i=0;i<nummat;i++)
	{
		tmp[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmp[i], 0, sizeof(doublecomplex)*size*size);
	}

	symmsc(mat_c22);
	for(i=rank;i<nummat;i+=numprocs)
	{
		symmsc2(mat_c2, i);
		for(j=0;j<=28;j+=4)
		{
			a=mat_c2[j];
			b=mat_c2[j+1];
			c=mat_c2[j+2];
			con=mat_c2[j+3];
			if(a>=0)
			{
				if(a==i)
				{
					for(k=0;k<=60;k+=4)
					{
						e=mat_c22[k];
						f=mat_c22[k+1];
						g=mat_c22[k+2];
						con2=mat_c22[k+3];
						if(e>=0)
						{
							if(b==e)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
								AntiComm3(tmp[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(f!=g)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
									AntiComm3(tmp[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(b==f)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
								AntiComm3(tmp[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
									AntiComm3(tmp[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(b==g)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
								AntiComm3(tmp[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(e!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
									AntiComm3(tmp[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							if(b!=c)
							{
								if(c==e)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
									AntiComm3(tmp[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(f!=g)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
										AntiComm3(tmp[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==f)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
									AntiComm3(tmp[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=e)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
										AntiComm3(tmp[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==g)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
									AntiComm3(tmp[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(e!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
										AntiComm3(tmp[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
							}
						}
						else
						{
							continue;
						}
					}
					Multi6(tmp[i], pmat[b], pmat[c], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					if(b!=c)
					{
						Multi6(tmp[i], pmat[c], pmat[b], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					}
				}
				else if(b==i)
				{
					for(k=0;k<=60;k+=4)
					{
						e=mat_c22[k];
						f=mat_c22[k+1];
						g=mat_c22[k+2];
						con2=mat_c22[k+3];

						if(e>=0)
						{
							if(a==e)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
								AntiComm3(tmp[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
									AntiComm3(tmp[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==f)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
								AntiComm3(tmp[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
									AntiComm3(tmp[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==g)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
								AntiComm3(tmp[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(f!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
									AntiComm3(tmp[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							if(a!=c)
							{
								if(c==e)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
									AntiComm3(tmp[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
										AntiComm3(tmp[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==f)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
									AntiComm3(tmp[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=e)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
										AntiComm3(tmp[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==g)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
									AntiComm3(tmp[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(f!=e)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
										AntiComm3(tmp[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
							}
						}
						else
						{
							continue;
						}
					}
					Multi6(tmp[i], pmat[a], pmat[c], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					if(a!=c)
					{
						Multi6(tmp[i], pmat[c], pmat[a], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					}
				}
				else if(c==i)
				{
					for(k=0;k<=60;k+=4)
					{
						e=mat_c22[k];
						f=mat_c22[k+1];
						g=mat_c22[k+2];
						con2=mat_c22[k+3];

						if(e>=0)
						{
							if(a==e)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
								AntiComm3(tmp[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
									AntiComm3(tmp[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==f)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
								AntiComm3(tmp[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
									AntiComm3(tmp[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==g)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
								AntiComm3(tmp[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(e!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
									AntiComm3(tmp[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							if(a!=b)
							{
								if(b==e)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
									AntiComm3(tmp[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
										AntiComm3(tmp[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(b==f)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
									AntiComm3(tmp[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(e!=g)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
										AntiComm3(tmp[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(b==g)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
									AntiComm3(tmp[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(e!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
										AntiComm3(tmp[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
							}
						}
						else
						{
							continue;
						}
					}
					Multi6(tmp[i], pmat[a], pmat[b], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					if(a!=b)
					{
						Multi6(tmp[i], pmat[b], pmat[a], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					}
				}
			}
			else
			{
				continue;
			}
		}
		AddMat7(tmp[i], pmat[i], -2.0*C3*C3*eps*lambda*alpha*alpha/sqrt(size), size);
	}

	for(i=0;i<nummat;i++)
	{
		MPI_Reduce(tmp[i], tmpmat[i], size*size, CMPI_COMPLEX_TYPE, CMPI_ADDCMAT_OP, 0, MPI_COMM_WORLD);
	}

	free(tmpmat1);
	free(tmpmat2);
	for(i=0;i<nummat;i++)
	{
		free(tmp[i]);
	}

	return 0;
}
