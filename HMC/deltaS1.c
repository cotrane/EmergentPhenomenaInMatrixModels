#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"

/*
 * variation of extra term in modified YM-Myers model; for more information see thesis, chapter 8
 */

int deltaS1(doublecomplex *pmat[], doublecomplex *tmpmat[], int nummat, int size, double eps, double C2, double C3, double lambda, double alpha)
{
	int i,j,k,a,b,c,e,f,g;
	double mat_c2[32], mat_c22[64];
	double con=0, con2=0;
	doublecomplex *tmpmat1, *tmpmat2;

	tmpmat1 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
	tmpmat2 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);

	// read in nonzero values of totally symmetric tensor d_{\mu\nu\rho}
	symmsc(mat_c22);
	for(i=0;i<nummat;i++)
	{
		// read in non-zero values of d_{\mu\nu\rho} for specific value of \mu
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
								AntiComm3(tmpmat[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(f!=g)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
									AntiComm3(tmpmat[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(b==f)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
								AntiComm3(tmpmat[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
									AntiComm3(tmpmat[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(b==g)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
								AntiComm3(tmpmat[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(e!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
									AntiComm3(tmpmat[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							if(b!=c)
							{
								if(c==e)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
									AntiComm3(tmpmat[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(f!=g)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
										AntiComm3(tmpmat[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==f)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
									AntiComm3(tmpmat[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=e)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
										AntiComm3(tmpmat[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==g)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
									AntiComm3(tmpmat[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(e!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
										AntiComm3(tmpmat[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
							}
						}
						else
						{
							continue;
						}
					}
					Multi6(tmpmat[i], pmat[b], pmat[c], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					if(b!=c)
					{
						Multi6(tmpmat[i], pmat[c], pmat[b], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
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
								AntiComm3(tmpmat[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
									AntiComm3(tmpmat[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==f)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
								AntiComm3(tmpmat[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
									AntiComm3(tmpmat[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==g)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
								AntiComm3(tmpmat[i], pmat[c], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(f!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
									AntiComm3(tmpmat[i], pmat[c], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							if(a!=c)
							{
								if(c==e)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
									AntiComm3(tmpmat[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
										AntiComm3(tmpmat[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==f)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
									AntiComm3(tmpmat[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=e)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
										AntiComm3(tmpmat[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(c==g)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
									AntiComm3(tmpmat[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(f!=e)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
										AntiComm3(tmpmat[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
							}
						}
						else
						{
							continue;
						}
					}
					Multi6(tmpmat[i], pmat[a], pmat[c], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					if(a!=c)
					{
						Multi6(tmpmat[i], pmat[c], pmat[a], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
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
								AntiComm3(tmpmat[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
									AntiComm3(tmpmat[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==f)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
								AntiComm3(tmpmat[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(g!=e)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
									AntiComm3(tmpmat[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							else if(a==g)
							{
								memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
								Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
								AntiComm3(tmpmat[i], pmat[b], tmpmat1, size, -2.0*con*con2*eps*lambda);
								if(e!=f)
								{
									memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
									AntiComm3(tmpmat[i], pmat[b], tmpmat2, size, -2.0*con*con2*eps*lambda);
								}
							}
							if(a!=b)
							{
								if(b==e)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[g], size, 1.0);
									AntiComm3(tmpmat[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(g!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[f], size, 1.0);
										AntiComm3(tmpmat[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(b==f)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[e], pmat[g], size, 1.0);
									AntiComm3(tmpmat[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(e!=g)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[g], pmat[e], size, 1.0);
										AntiComm3(tmpmat[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
								else if(b==g)
								{
									memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
									Multi5(tmpmat1, pmat[f], pmat[e], size, 1.0);
									AntiComm3(tmpmat[i], pmat[a], tmpmat1, size, -2.0*con*con2*eps*lambda);
									if(e!=f)
									{
										memset(tmpmat2, 0, sizeof(doublecomplex)*size*size);
										Multi5(tmpmat2, pmat[e], pmat[f], size, 1.0);
										AntiComm3(tmpmat[i], pmat[a], tmpmat2, size, -2.0*con*con2*eps*lambda);
									}
								}
							}
						}
						else
						{
							continue;
						}
					}
					Multi6(tmpmat[i], pmat[a], pmat[b], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					if(a!=b)
					{
						Multi6(tmpmat[i], pmat[b], pmat[a], size, 6.0*C3*con*alpha*eps*lambda/pow(size,0.25));
					}
				}
			}
			else
			{
				continue;
			}
		}
		AddMat7(tmpmat[i], pmat[i], -2.0*C3*C3*eps*lambda*alpha*alpha/sqrt(size), size);
	}

	free(tmpmat1);
	free(tmpmat2);

	return 0;
}
