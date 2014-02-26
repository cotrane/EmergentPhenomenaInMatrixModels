#include<stdlib.h>
#include<f2c.h>
#include<stdio.h>
#include<math.h>
#include "hmc.h"
#include "montecarlo.h"
#include "RandomGens.h"

/*
 * generate gaussian, diagonal, traceless matrix pmom; parameter 'a' allows to increase the initial mean momentum
 */
int gen_gaussmomcplx(double pmom[], int MATRIX_SIZE, double a)
{
	int m;
	double lastel=0.0;

	for(m = 0; m < (MATRIX_SIZE-1) ;m++)
	{
		pmom[m] = a*gauss_randomnr(1);
		lastel += pmom[m];
	}
	pmom[(MATRIX_SIZE-1)] = -lastel;

	return 0;
}

/*
 * generates gaussian, diagonal matrix pmom; parameter 'a' allows to change initial mean momentum
 */
int gen_gaussmomcplxtrace(double pmom[], int MATRIX_SIZE, double a)
{
	int m;

	for(m = 0; m < (MATRIX_SIZE) ;m++)
	{
		pmom[m] = a*gauss_randomnr(1);
	}

	return 0;
}
