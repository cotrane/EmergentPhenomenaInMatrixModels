#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<time.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include "montecarlo.h"

int symmsc(double *mat_c2)
{
	mat_c2[0]=0;
	mat_c2[1]=0;
	mat_c2[2]=7;
	mat_c2[3]=1.0/sqrt(3.0);

	mat_c2[4]=0;
	mat_c2[5]=3;
	mat_c2[6]=5;
	mat_c2[7]=0.5;

	mat_c2[8]=0;
	mat_c2[9]=4;
	mat_c2[10]=6;
	mat_c2[11]=0.5;

	mat_c2[12]=1;
	mat_c2[13]=1;
	mat_c2[14]=7;
	mat_c2[15]=1.0/sqrt(3.0);

	mat_c2[16]=1;
	mat_c2[17]=3;
	mat_c2[18]=6;
	mat_c2[19]=-0.5;

	mat_c2[20]=1;
	mat_c2[21]=4;
	mat_c2[22]=5;
	mat_c2[23]=0.5;

	mat_c2[24]=2;
	mat_c2[25]=2;
	mat_c2[26]=7;
	mat_c2[27]=1.0/sqrt(3.0);

	mat_c2[28]=2;
	mat_c2[29]=3;
	mat_c2[30]=3;
	mat_c2[31]=0.5;;

	mat_c2[32]=2;
	mat_c2[33]=4;
	mat_c2[34]=4;
	mat_c2[35]=0.5;

	mat_c2[36]=2;
	mat_c2[37]=5;
	mat_c2[38]=5;
	mat_c2[39]=-0.5;

	mat_c2[40]=2;
	mat_c2[41]=6;
	mat_c2[42]=6;
	mat_c2[43]=-0.5;

	mat_c2[44]=3;
	mat_c2[45]=3;
	mat_c2[46]=7;
	mat_c2[47]=-1.0/(2.0*sqrt(3.0));

	mat_c2[48]=4;
	mat_c2[49]=4;
	mat_c2[50]=7;
	mat_c2[51]=-1.0/(2.0*sqrt(3.0));

	mat_c2[52]=5;
	mat_c2[53]=5;
	mat_c2[54]=7;
	mat_c2[55]=-1.0/(2.0*sqrt(3.0));

	mat_c2[56]=6;
	mat_c2[57]=6;
	mat_c2[58]=7;
	mat_c2[59]=-1.0/(2.0*sqrt(3.0));

	mat_c2[60]=7;
	mat_c2[61]=7;
	mat_c2[62]=7;
	mat_c2[63]=-1.0/(sqrt(3.0));

	return 0;
}

int symmsc2(double *mat_c2, int mat)
{
	if(mat==0)				//labels different matrices in the array
	{
		mat_c2[0]=0;
		mat_c2[1]=0;
		mat_c2[2]=7;
		mat_c2[3]=1.0/sqrt(3.0);

		mat_c2[4]=0;
		mat_c2[5]=3;
		mat_c2[6]=5;
		mat_c2[7]=0.5;

		mat_c2[8]=0;
		mat_c2[9]=4;
		mat_c2[10]=6;
		mat_c2[11]=0.5;

		mat_c2[12]=-1;
		mat_c2[13]=-1;
		mat_c2[14]=-1;
		mat_c2[15]=-1;

		mat_c2[16]=-1;
		mat_c2[17]=-1;
		mat_c2[18]=-1;
		mat_c2[19]=-1;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==1)
	{
		mat_c2[0]=1;
		mat_c2[1]=1;
		mat_c2[2]=7;
		mat_c2[3]=1.0/sqrt(3.0);

		mat_c2[4]=1;
		mat_c2[5]=3;
		mat_c2[6]=6;
		mat_c2[7]=-0.5;

		mat_c2[8]=1;
		mat_c2[9]=4;
		mat_c2[10]=5;
		mat_c2[11]=0.5;

		mat_c2[12]=-1;
		mat_c2[13]=-1;
		mat_c2[14]=-1;
		mat_c2[15]=-1;

		mat_c2[16]=-1;
		mat_c2[17]=-1;
		mat_c2[18]=-1;
		mat_c2[19]=-1;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==2)
	{
		mat_c2[0]=2;
		mat_c2[1]=2;
		mat_c2[2]=7;
		mat_c2[3]=1.0/sqrt(3);

		mat_c2[4]=2;
		mat_c2[5]=3;
		mat_c2[6]=3;
		mat_c2[7]=0.5;;

		mat_c2[8]=2;
		mat_c2[9]=4;
		mat_c2[10]=4;
		mat_c2[11]=0.5;

		mat_c2[12]=2;
		mat_c2[13]=5;
		mat_c2[14]=5;
		mat_c2[15]=-0.5;

		mat_c2[16]=2;
		mat_c2[17]=6;
		mat_c2[18]=6;
		mat_c2[19]=-0.5;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==3)
	{
		mat_c2[0]=3;
		mat_c2[1]=3;
		mat_c2[2]=7;
		mat_c2[3]=-1.0/(2.0*sqrt(3.0));

		mat_c2[4]=1;
		mat_c2[5]=3;
		mat_c2[6]=6;
		mat_c2[7]=-0.5;

		mat_c2[8]=2;
		mat_c2[9]=3;
		mat_c2[10]=3;
		mat_c2[11]=0.5;

		mat_c2[12]=0;
		mat_c2[13]=3;
		mat_c2[14]=5;
		mat_c2[15]=0.5;

		mat_c2[16]=-1;
		mat_c2[17]=-1;
		mat_c2[18]=-1;
		mat_c2[19]=-1;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==4)
	{
		mat_c2[0]=4;
		mat_c2[1]=4;
		mat_c2[2]=7;
		mat_c2[3]=-1.0/(2.0*sqrt(3.0));

		mat_c2[4]=0;
		mat_c2[5]=4;
		mat_c2[6]=6;
		mat_c2[7]=0.5;

		mat_c2[8]=1;
		mat_c2[9]=4;
		mat_c2[10]=5;
		mat_c2[11]=0.5;

		mat_c2[12]=2;
		mat_c2[13]=4;
		mat_c2[14]=4;
		mat_c2[15]=0.5;

		mat_c2[16]=-1;
		mat_c2[17]=-1;
		mat_c2[18]=-1;
		mat_c2[19]=-1;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==5)
	{
		mat_c2[0]=5;
		mat_c2[1]=5;
		mat_c2[2]=7;
		mat_c2[3]=-1.0/(2.0*sqrt(3.0));

		mat_c2[4]=0;
		mat_c2[5]=3;
		mat_c2[6]=5;
		mat_c2[7]=0.5;

		mat_c2[8]=1;
		mat_c2[9]=4;
		mat_c2[10]=5;
		mat_c2[11]=0.5;

		mat_c2[12]=2;
		mat_c2[13]=5;
		mat_c2[14]=5;
		mat_c2[15]=-0.5;

		mat_c2[16]=-1;
		mat_c2[17]=-1;
		mat_c2[18]=-1;
		mat_c2[19]=-1;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==6)
	{
		mat_c2[0]=6;
		mat_c2[1]=6;
		mat_c2[2]=7;
		mat_c2[3]=-1.0/(2.0*sqrt(3.0));

		mat_c2[4]=0;
		mat_c2[5]=4;
		mat_c2[6]=6;
		mat_c2[7]=0.5;

		mat_c2[8]=1;
		mat_c2[9]=3;
		mat_c2[10]=6;
		mat_c2[11]=-0.5;

		mat_c2[12]=2;
		mat_c2[13]=6;
		mat_c2[14]=6;
		mat_c2[15]=-0.5;

		mat_c2[16]=-1;
		mat_c2[17]=-1;
		mat_c2[18]=-1;
		mat_c2[19]=-1;

		mat_c2[20]=-1;
		mat_c2[21]=-1;
		mat_c2[22]=-1;
		mat_c2[23]=-1;

		mat_c2[24]=-1;
		mat_c2[25]=-1;
		mat_c2[26]=-1;
		mat_c2[27]=-1;

		mat_c2[28]=-1;
		mat_c2[29]=-1;
		mat_c2[30]=-1;
		mat_c2[31]=-1;
	}
	else if(mat==7)
	{
		mat_c2[0]=0;
		mat_c2[1]=0;
		mat_c2[2]=7;
		mat_c2[3]=1.0/(sqrt(3.0));

		mat_c2[4]=1;
		mat_c2[5]=1;
		mat_c2[6]=7;
		mat_c2[7]=1.0/(sqrt(3.0));

		mat_c2[8]=2;
		mat_c2[9]=2;
		mat_c2[10]=7;
		mat_c2[11]=1.0/(sqrt(3.0));

		mat_c2[12]=7;
		mat_c2[13]=7;
		mat_c2[14]=7;
		mat_c2[15]=-1.0/(sqrt(3.0));

		mat_c2[16]=3;
		mat_c2[17]=3;
		mat_c2[18]=7;
		mat_c2[19]=-1.0/(2.0*sqrt(3.0));

		mat_c2[20]=4;
		mat_c2[21]=4;
		mat_c2[22]=7;
		mat_c2[23]=-1.0/(2.0*sqrt(3.0));

		mat_c2[24]=5;
		mat_c2[25]=5;
		mat_c2[26]=7;
		mat_c2[27]=-1.0/(2.0*sqrt(3.0));

		mat_c2[28]=6;
		mat_c2[29]=6;
		mat_c2[30]=7;
		mat_c2[31]=-1.0/(2.0*sqrt(3.0));
	}

	return 0;
}

