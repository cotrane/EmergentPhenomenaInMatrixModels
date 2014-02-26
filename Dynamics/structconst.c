#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"

/*
 * structure constant f_{abc} for 8MM
 */
int structconst(double *mat_c)
{
	mat_c[0]=0;
	mat_c[1]=1;
	mat_c[2]=2;
	mat_c[3]=1.0;

	mat_c[4]=0;
	mat_c[5]=3;
	mat_c[6]=6;
	mat_c[7]=0.5;

	mat_c[8]=0;
	mat_c[9]=4;
	mat_c[10]=5;
	mat_c[11]=-0.5;

	mat_c[12]=1;
	mat_c[13]=3;
	mat_c[14]=5;
	mat_c[15]=0.5;

	mat_c[16]=1;
	mat_c[17]=4;
	mat_c[18]=6;
	mat_c[19]=0.5;

	mat_c[20]=2;
	mat_c[21]=3;
	mat_c[22]=4;
	mat_c[23]=0.5;

	mat_c[24]=2;
	mat_c[25]=5;
	mat_c[26]=6;
	mat_c[27]=-0.5;

	mat_c[28]=3;
	mat_c[29]=4;
	mat_c[30]=7;
	mat_c[31]=sqrt(3.0)/2.0;

	mat_c[32]=5;
	mat_c[33]=6;
	mat_c[34]=7;
	mat_c[35]=sqrt(3.0)/2.0;

	return 0;
}

/*
 * thecode.c
 *
 *  Created on: 29 Sep 2010
 *      Author: tkaltenbrunner
 */

