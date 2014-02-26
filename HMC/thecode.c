#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"

/*
 * definition of structure constant f_{\mu\nu\rho} for 8MM to read in specific non-zero values for certain matrix number 'mat'
 */
int thecode(int mat, double *mat_c)
{
					if(mat==0)			//labels different matrices in the array
					{
						mat_c[0]=1;			//value for second entry in commutator
						mat_c[1]=2;			//value for first entry in commutator
						mat_c[2]=1.0;
						mat_c[3]=0;
						mat_c[4]=2;
						mat_c[5]=1;

						mat_c[6]=3;
						mat_c[7]=6;
						mat_c[8]=0.5;
						mat_c[9]=3;
						mat_c[10]=5;
						mat_c[11]=4;

						mat_c[12]=4;
						mat_c[13]=5;
						mat_c[14]=-0.5;
						mat_c[15]=6;
						mat_c[16]=8;
						mat_c[17]=7;

						mat_c[18]=-1;
						mat_c[19]=-1;
						mat_c[20]=-1;
						mat_c[21]=-1;
						mat_c[22]=-1;
						mat_c[23]=-1;
					}
					else if(mat==1)
					{
						mat_c[0]=2;			//value for second entry in commutator
						mat_c[1]=0;			//value for first entry in commutator
						mat_c[2]=1.0;
						mat_c[3]=1;
						mat_c[4]=0;
						mat_c[5]=2;

						mat_c[6]=3;
						mat_c[7]=5;
						mat_c[8]=0.5;
						mat_c[9]=9;
						mat_c[10]=11;
						mat_c[11]=10;

						mat_c[12]=4;
						mat_c[13]=6;
						mat_c[14]=0.5;
						mat_c[15]=12;
						mat_c[16]=14;
						mat_c[17]=13;

						mat_c[18]=-1;
						mat_c[19]=-1;
						mat_c[20]=-1;
						mat_c[21]=-1;
						mat_c[22]=-1;
						mat_c[23]=-1;
					}
					else if(mat==2)
					{
						mat_c[0]=0;			//value for second entry in commutator
						mat_c[1]=1;			//value for first entry in commutator
						mat_c[2]=1.0;
						mat_c[3]=2;
						mat_c[4]=1;
						mat_c[5]=0;

						mat_c[6]=3;
						mat_c[7]=4;
						mat_c[8]=0.5;
						mat_c[9]=15;
						mat_c[10]=17;
						mat_c[11]=16;

						mat_c[12]=5;
						mat_c[13]=6;
						mat_c[14]=-0.5;
						mat_c[15]=18;
						mat_c[16]=20;
						mat_c[17]=19;

						mat_c[18]=-1;
						mat_c[19]=-1;
						mat_c[20]=-1;
						mat_c[21]=-1;
						mat_c[22]=-1;
						mat_c[23]=-1;
					}
					else if(mat==3)
					{
						mat_c[0]=6;			//value for second entry in commutator
						mat_c[1]=0;			//value for first entry in commutator
						mat_c[2]=0.5;
						mat_c[3]=4;
						mat_c[4]=3;
						mat_c[5]=5;

						mat_c[6]=5;
						mat_c[7]=1;
						mat_c[8]=0.5;
						mat_c[9]=10;
						mat_c[10]=9;
						mat_c[11]=11;

						mat_c[12]=4;
						mat_c[13]=2;
						mat_c[14]=0.5;
						mat_c[15]=16;
						mat_c[16]=15;
						mat_c[17]=17;

						mat_c[18]=4;
						mat_c[19]=7;
						mat_c[20]=sqrt(3.0)/2.0;
						mat_c[21]=21;
						mat_c[22]=23;
						mat_c[23]=22;
					}
					else if(mat==4)
					{
						mat_c[0]=5;			//value for second entry in commutator
						mat_c[1]=0;			//value for first entry in commutator
						mat_c[2]=-0.5;
						mat_c[3]=7;
						mat_c[4]=6;
						mat_c[5]=8;

						mat_c[6]=6;
						mat_c[7]=1;
						mat_c[8]=0.5;
						mat_c[9]=13;
						mat_c[10]=12;
						mat_c[11]=14;

						mat_c[12]=2;
						mat_c[13]=3;
						mat_c[14]=0.5;
						mat_c[15]=17;
						mat_c[16]=16;
						mat_c[17]=15;

						mat_c[18]=7;
						mat_c[19]=3;
						mat_c[20]=sqrt(3.0)/2.0;
						mat_c[21]=22;
						mat_c[22]=21;
						mat_c[23]=23;
					}
					else if(mat==5)
					{
						mat_c[0]=0;			//value for second entry in commutator
						mat_c[1]=4;			//value for first entry in commutator
						mat_c[2]=-0.5;
						mat_c[3]=8;
						mat_c[4]=7;
						mat_c[5]=6;

						mat_c[6]=1;
						mat_c[7]=3;
						mat_c[8]=0.5;
						mat_c[9]=11;
						mat_c[10]=10;
						mat_c[11]=9;

						mat_c[12]=6;
						mat_c[13]=2;
						mat_c[14]=-0.5;
						mat_c[15]=19;
						mat_c[16]=18;
						mat_c[17]=20;

						mat_c[18]=6;
						mat_c[19]=7;
						mat_c[20]=sqrt(3.0)/2.0;
						mat_c[21]=24;
						mat_c[22]=26;
						mat_c[23]=25;
					}
					else if(mat==6)
					{
						mat_c[0]=0;			//value for second entry in commutator
						mat_c[1]=3;			//value for first entry in commutator
						mat_c[2]=0.5;
						mat_c[3]=5;
						mat_c[4]=4;
						mat_c[5]=3;

						mat_c[6]=1;
						mat_c[7]=4;
						mat_c[8]=0.5;
						mat_c[9]=14;
						mat_c[10]=13;
						mat_c[11]=12;

						mat_c[12]=2;
						mat_c[13]=5;
						mat_c[14]=-0.5;
						mat_c[15]=20;
						mat_c[16]=19;
						mat_c[17]=18;

						mat_c[18]=7;
						mat_c[19]=5;
						mat_c[20]=sqrt(3.0)/2.0;
						mat_c[21]=25;
						mat_c[22]=24;
						mat_c[23]=26;
					}
					else if(mat==7)
					{
						mat_c[0]=3;			//value for second entry in commutator
						mat_c[1]=4;			//value for first entry in commutator
						mat_c[2]=sqrt(3.0)/2.0;
						mat_c[3]=23;
						mat_c[4]=22;
						mat_c[5]=21;

						mat_c[6]=5;
						mat_c[7]=6;
						mat_c[8]=sqrt(3.0)/2.0;
						mat_c[9]=26;
						mat_c[10]=25;
						mat_c[11]=24;

						mat_c[12]=-1;
						mat_c[13]=-1;
						mat_c[14]=-1;
						mat_c[15]=-1;
						mat_c[16]=-1;
						mat_c[17]=-1;

						mat_c[18]=-1;
						mat_c[19]=-1;
						mat_c[20]=-1;
						mat_c[21]=-1;
						mat_c[22]=-1;
						mat_c[23]=-1;
					}

					return 0;
}

/*
 * definition of structure constant f_{\mu\nu\rho} for 8MM
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

