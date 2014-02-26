#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"

/*
 * computes value of action S[mass, pmat] for matrices with trace
 */
double actionYM(double pmat[], int size, double mass)
{
	int i,j;
	double trYM=0;

	//YM-Term

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			if(j!=i)
			{
				trYM += -0.5*log((pmat[i]-pmat[j])*(pmat[i]-pmat[j])) + 0.5*log(size+size*(pmat[i]-pmat[j])*(pmat[i]-pmat[j])/(2.*mass*mass));
			}
		}
	}

	//Mass-Term

	for(i=0;i<size;i++)
	{
		trYM += size*pmat[i]*pmat[i];
	}

	return trYM;
}

/*
 * computes value of action S[g,pmat] for matrices with trace
 */
double actionYM_N(double pmat[], int size, double g)
{
	int i,j;
	double trYM=0;

	//YM-Term

	if(g!=0){
		for(i=0;i<size;i++)
		{
			for(j=0;j<size;j++)
			{
				if(j!=i)
				{
					trYM += -0.5*log((pmat[i]-pmat[j])*(pmat[i]-pmat[j])) + 0.5*log(1+g*g*(pmat[i]-pmat[j])*(pmat[i]-pmat[j]));
				}
			}
		}
	}

	//Mass-Term

	for(i=0;i<size;i++)
	{
		trYM += pmat[i]*pmat[i];
	}

	return trYM;
}

