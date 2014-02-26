/*
 * nonHermitianMatrix.c
 *
 *  Created on: 23 Jan 2013
 *      Author: tkaltenbrunner
 */
#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<math.h>
#include"hmc.h"
#include "RandomGens.h"

int genNonHermitianMatrix(doublecomplex *pmat, int size)
{
	int m, n;
	doublecomplex lastel;

	lastel.r=0; lastel.i=0;
	for(m = 0; m < (size-1) ;m++)
	{
		for(n=m; n < size ; n++)
		{
			if(m == n)
			{
				pmat[m*size + n].r = 1./(sqrt(2.0)) * gauss_randomnr(1);//prng_get_double();
				pmat[m*size + n].i = 1./(sqrt(2.0)) * gauss_randomnr(2);//prng_get_double();
				lastel.r += pmat[m*size + n].r;
				lastel.i += pmat[m*size + n].i;
			}
			else
			{
				pmat[m*size + n].r = 1./(sqrt(2.0)) * gauss_randomnr(1);//prng_get_double();
				pmat[m*size + n].i = 1./(sqrt(2.0)) * gauss_randomnr(2);//prng_get_double();
				pmat[n*size + m].r = 1./(sqrt(2.0)) * gauss_randomnr(1);//prng_get_double();
				pmat[n*size + m].i = 1./(sqrt(2.0)) * gauss_randomnr(2);//prng_get_double();
			}
		}
	}
	pmat[(size-1)*size + (size-1)].r = -lastel.r;
	pmat[(size-1)*size + (size-1)].i = -lastel.i;

	return 0;
}
