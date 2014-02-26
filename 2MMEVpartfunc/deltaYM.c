#include<stdlib.h>
#include<stdio.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "hmc.h"

/*
 * adds \frac{\delta S[mass, pmat]}{\delta pmat} to pmom
 */
int deltaYM(double pmat[], double pmom[], int size, double eps, double mass)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			if(j!=i)
			{
				pmom[i] += 2.*eps*( 1. / ((pmat[i]-pmat[j])*(1.+(pmat[i]-pmat[j])*(pmat[i]-pmat[j])/(2.*mass*mass))));
			}
		}
		pmom[i] += -2.*eps*size*pmat[i];
	}

	return 0;
}

/*
 * adds \frac{\delta S[g, pmat]}{\delta pmat} to pmom
 */
int deltaYM_N(double pmat[], double pmom[], int size, double eps, double g)
{
	int i,j;

	for(i=0;i<size;i++)
	{
		if(g!=0){
			for(j=0;j<size;j++)
			{
				if(j!=i)
				{
					pmom[i] += 2.*eps*( 1. / ((pmat[i]-pmat[j])*(1+g*g*(pmat[i]-pmat[j])*(pmat[i]-pmat[j]))));
				}
			}
		}
		pmom[i] += -2.*eps*pmat[i];
	}

	return 0;
}
