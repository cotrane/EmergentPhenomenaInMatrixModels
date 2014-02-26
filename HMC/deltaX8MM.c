#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "prng.h"
#include "montecarlo.h"
#include "hmc.h"
#include "MatrixMan.h"

/*
 * variation of action for 8MM and modified 8MM which is added to momentum matrices pmom
 */

int deltaX8MM(doublecomplex *pmat[], doublecomplex *pmom[], int nummat, int size, double coupling, double eps, double C2, double C3, double lambda)
{
	int i,j,a,b,c,save, savea=0, saveb=0, savec=0, pos1, pos2, turn=0;
	double mat_c[36], con, tr=0;
	doublecomplex *pcomm[28], *tmpmat[nummat];

	for(i=0;i<nummat;i++)
	{
		tmpmat[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(tmpmat[i], 0, sizeof(doublecomplex)*size*size);
	}
	for(i=0;i<28;i++)
	{
		pcomm[i] = (doublecomplex*) malloc(sizeof(doublecomplex)*size*size);
		memset(pcomm[i], 0, sizeof(doublecomplex)*size*size);
	}

	// if lambda not zero computation of variation of additional term in modified 8MM (see thesis chapter 8)
	if(lambda!=0.0)
	{
		deltaS1(pmat, tmpmat, nummat, size, eps, C2, C3, lambda, coupling);
	}

	//computes commutator
	for(a=0; a<(nummat-1); a++)
	{
		pos1=0;
		for(i=0; i<a; i++)
		{
			pos1+=((nummat-1)-i);
		}
		pos2=0;
		for(b=(a+1); b<nummat; b++)
		{
			Comm(pcomm[pos1+pos2], pmat[a], pmat[b], size, 1.0);
			pos2++;
		}
	}
	// computes YM Term
	for(a=0;a<nummat;a++)
	{
		for(b=0;b<nummat;b++)
		{
			if(b!=a)
			{
				if(a<b)
				{
					con=1.0;
					pos1=0;
					for(i=0; i<a; i++)
					{
						pos1+=((nummat-1)-i);
					}
					Commend2(tmpmat[a], pmat[b], pcomm[pos1+b-1-a], size, eps*size*con);
				}
				if(b<a)
				{
					con=-1.0;
					pos1=0;
					for(i=0; i<b; i++)
					{
						pos1+=((nummat-1)-i);
					}
					Commend2(tmpmat[a], pmat[b], pcomm[pos1+a-1-b], size, eps*size*con);
				}

			}
		}
	}

	//Myers-Term
	structconst(mat_c);
	for(i=0;i<=32;i+=4)
	{
		a=mat_c[i];
		b=mat_c[i+1];
		c=mat_c[i+2];
		con=mat_c[i+3];

		savea = a;
		//abc
		pos1=0;
		for(j=0;j<b;j++)
		{
			pos1+=((nummat-1)-j);
		}
		AddMat4Dyn(tmpmat[a], pcomm[pos1+c-1-b], size, -2.0*eps*coupling*con*size/pow(size,0.25));

		//bca
		save = a;
		a=b;
		b=c;
		c=save;

		if(b>c)
		{
			con = -con;
			saveb=b;
			savec=c;
			b=c;
			c=saveb;
			turn=1;
		}

		pos1=0;
		for(j=0;j<b;j++)
		{
			pos1+=((nummat-1)-j);
		}
		AddMat4Dyn(tmpmat[a], pcomm[pos1+c-1-b], size, -2.0*eps*coupling*con*size/pow(size,0.25));

		if(turn==1)
		{
			b=saveb;
			c=savec;
			con=-con;
			turn=0;
		}

		//cab
		save = a;
		a=b;
		b=c;
		c=save;

		if(b>c)
		{
			con = -con;
			saveb=b;
			savec=c;
			b=c;
			c=saveb;
			turn=1;
		}

		pos1=0;
		for(j=0;j<b;j++)
		{
			pos1+=((nummat-1)-j);
		}
		AddMat4Dyn(tmpmat[a], pcomm[pos1+c-1-b], size, -2.0*eps*coupling*con*size/pow(size,0.25));
	}

	// implementation of tracelessness condition where X_{NN} = - \sum_{i=0}^{N-1} X_{ii}
	for(j=0;j<nummat;j++)
	{
		tr=0;
		for(i=0;i<(size-1);i++)
		{
			(tmpmat[j]+i*size+i)->r -= (tmpmat[j]+(size-1)*size+(size-1))->r;
			tr += (tmpmat[j]+i*size+i)->r;
		}
		(tmpmat[j]+(size-1)*size+(size-1))->r = -tr;
	}

	// adding results to momentum matrices pmom
	for(i=0;i<nummat;i++)
	{
		AddMat6(pmom[i], tmpmat[i], 1.0, size);
	}

	for(i=0;i<nummat;i++)
	{
		free(tmpmat[i]);
	}
	for(i=0;i<28;i++)
	{
		free(pcomm[i]);
	}

	return 0;
}
