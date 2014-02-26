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
 * compute contribution to action from additional term in modified 8MM (see chapter 8, thesis)
 */
double S1(doublecomplex *pmat[], int size, double C2, double C3, double lambda, double alpha, int nummat)
{
	int a,b,c,i,j;
	double con, trrand=0;
	double mat_c2[32];
	doublecomplex *tmpmat1;

	tmpmat1 = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);

	for(j=0;j<nummat;j++)
	{
		memset(tmpmat1, 0, sizeof(doublecomplex)*size*size);
		// read in non-zero values of structure constant that include specific matrix number 'j'
		symmsc2(mat_c2, j);
		for(i=0;i<=28;i+=4)
		{
			a=mat_c2[i];
			b=mat_c2[i+1];
			c=mat_c2[i+2];
			con=mat_c2[i+3];
			if(a>=0)
			{
				if(a==j)
				{
					if(b!=c)
					{
						AntiComm4(tmpmat1, pmat[b], pmat[c], size, con);
					}
					else
					{
						Multi5(tmpmat1, pmat[b], pmat[b], size, con);
					}
				}
				else if(b==j)
				{
					if(a!=c)
					{
						AntiComm4(tmpmat1, pmat[a], pmat[c], size, con);
					}
					else
					{
						Multi5(tmpmat1, pmat[a], pmat[a], size, con);
					}
				}
				else if(c==j)
				{
					if(b!=a)
					{
						AntiComm4(tmpmat1, pmat[b], pmat[a], size, con);
					}
					else
					{
						Multi5(tmpmat1, pmat[a], pmat[a], size, con);
					}
				}
			}
			else
			{
				continue;
			}
		}
		AddMat6(tmpmat1, pmat[j], -C3*alpha/pow(size, 0.25), size);
		diagMulti(&trrand, tmpmat1, tmpmat1, size, lambda);
	}

	free(tmpmat1);

	return trrand;
}
