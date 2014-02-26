#include<stdlib.h>
#include<math.h>
#include <f2c.h>
#include<clapack.h>
#include "mtwist.h"

	/*
	 * generate the gaussian distributed variables by Box-Muller-procedure:
	 * generate normal distribution out of flat distribution between [0,1)
	 * are distributed by P(xi) = 1/sqrt(2pi) exp(-xi^2/2)
	 * random eigenvalues between [0,1) are generated using Mersenne twistor
	 */


double gauss_randomnr(int n)
{
	double r, phi, grand=0;

	phi = 2.0*M_PI*(1.0-mt_ldrand());
	r = sqrt(-2.0*log(1.0-mt_ldrand()));
	if(n==1)
	{
		grand = cos(phi) *r;
	}
	if(n==2)
	{
		grand = sin(phi)*r;
	}
	return grand;
}

double gauss_randvarmean(int n, double var, double mean)
{
	double pii = 3.14159265;
	double r, phi, grand=0;

	phi = 2.0*pii*(1.0-mt_ldrand());
	r = sqrt(-2.0*log(1.0-mt_ldrand()));
	if(n==1)
	{
		grand = cos(phi)*r * var + mean;
	}
	if(n==2)
	{
		grand = sin(phi)*r * var + mean;
	}
	return grand;
}

doublecomplex gauss_randcomplex()
{
	double r, phi;
	doublecomplex grand;

	phi = 2.0*M_PI*(1.0-mt_ldrand());
	r = sqrt(-2.0*log(1.0-mt_ldrand()));

	grand.r = cos(phi)*r ;
	grand.i = sin(phi)*r;

	return grand;
}
