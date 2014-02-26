/*
 * GenGaussMatrix.c
 *
 *  Created on: 12 Jun 2012
 *      Author: tkaltenbrunner
 */
#include<stdio.h>
#include<stdlib.h>
#include"hmc.h"
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"
#include "MatrixMan.h"
#include "RandomGens.h"
#include "EV.h"

// generate matrix pmatXYZ[1] as gaussian matrix using matrix pmatXYZ[0] (see thesis)
int GenGaussMatrixY(doublecomplex *pout, double *pmat, double mass, int size)
{
	int m, n;
	double a=0, lastel=0;
	doublecomplex rannum;

	for(m = 0; m < (size) ;m++)
	{
		for(n=m; n < size ; n++)
		{
			if(m == n)
			{
				pout[m*size + n].r = 1/sqrt(2.0*size) * gauss_randcomplex().r;
				pout[m*size + n].i = 0.0;
			}
			else
			{
				a = 1.0/sqrt(2.0*size + 2.0*size*(pmat[m]-pmat[n])*(pmat[m]-pmat[n])/(2.0*mass*mass));
				rannum = gauss_randcomplex();
				pout[m*size + n].r = a*rannum.r/sqrt(2.0);
				pout[m*size + n].i = a*rannum.i/sqrt(2.0);

				pout[n*size + m].r = pout[m*size + n].r;
				pout[n*size + m].i = -pout[m*size + n].i;
			}
		}
	}

	return 0;
}
