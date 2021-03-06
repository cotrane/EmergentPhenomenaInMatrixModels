#include<stdio.h>
#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include<string.h>
#include<math.h>
#include"MatrixMan.h"
#include"RepsSUn.h"


int ScanMatrixZero(doublecomplex *pmat, int size);
int ScanMatrixDiag(doublecomplex *pmat, int size);

/*
 * code generates \gamma-matrices corresponding to dimensionality of the model defined by nummat
 * matrices are stored in pgamma[]
 */

int gammamatrices(doublecomplex *pgamma[], int nummat)
{
	int i,j,k;
	int N=0,n, ntmp,d, matsize=0, check=0;
	doublecomplex *gens[3], **tmp, **tmp2, *tmp3, *test;

	tmp = (doublecomplex**) malloc(sizeof(doublecomplex*)*3);
	tmp2 = (doublecomplex**) malloc(sizeof(doublecomplex*)*3);
	for(j=0;j<3;j++)
	{
		tmp[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*4);
		memset(tmp[j], 0, sizeof(doublecomplex)*4);
		tmp2[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*4);
		memset(tmp2[j], 0, sizeof(doublecomplex)*4);
		gens[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*4);
		memset(gens[j], 0, sizeof(doublecomplex)*4);
	}
	tmp3 = (doublecomplex*) malloc(sizeof(doublecomplex)*4);
	memset(tmp3, 0, sizeof(doublecomplex)*4);

	RepsSUn(gens, 2, 2, 3, 2.0);

	/*
	 * if there are only 2 matrices gammamatrices correspond to first 2 pauli matrices, generated by RepsSUn; in this case only
	 * next if-loop is executed
	 */
	if(nummat==2)
	{
		for(i=0;i<2;i++)
		{
			memcpy(pgamma[i], gens[i], sizeof(doublecomplex)*4);
		}
		return 0;
	}

	//define number of gamma matrices; 2*N+n - dimensional
	for(i=1;i<20;i++)
	{
		for(j=1;j<4;j++)
		{
			if((nummat/(2.*i+j)) == 1.0)
			{
				N=i;
				n=j;
				break;
			}
		}
		if(N!=0)
			break;
	}

	/*
	 * generating \gamma-matrices step by step; at each iteration dimension is increased by one
	 */
	for(i=1;i<(N+1);i++)
	{
		matsize=pow(2,i+1);
		d=2*i+3;
		if(d<(2*N+n))
			ntmp=2;
		else
			ntmp=n-1;

		for(j=0;j<(2*(i-1)+3);j++)
		{
			free(tmp2[j]);
		}
		tmp2 = (doublecomplex**) malloc(sizeof(doublecomplex*)*(2*i+(ntmp+1)));
		for(j=0;j<(2*i+ntmp+1);j++)
		{
			tmp2[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*matsize*matsize);
			memset(tmp2[j], 0, sizeof(doublecomplex)*matsize*matsize);
		}

		//special case where i=1 -> for n=0 no cross product; for n!=0 empty tmp-matrices!!
		if(i==1)
		{
			if(!ntmp)
			{
				for(j=0;j<2;j++)
					memcpy(tmp[j], gens[j], sizeof(doublecomplex)*4);
				AddMat8(tmp[2], gens[2], -1.0, 2);
				break;
			}
			else
			{
				//cross product with d-ntmp gammas and sigma1
				for(k=0;k<(2*i+1);k++)
				{
					Multilambda(tmp2[k], gens[k], gens[0], 2, 2);
				}
				//cross product with unity and sigma2 (and sigma 3 if ntmp=2)
				for(k=0;k<2;k++)
				{
					tmp3[k*(int)(pow(2,i))+k].r = 1.0;
				}
				for(j=(2*i+1);j<(2*i+1+ntmp);j++)
				{
					if(j==(d-2))
						Multilambda(tmp2[j], tmp3, gens[1], 2, 2);
					if(j==(d-1))
						Multilambda(tmp2[j], tmp3, gens[2], 2, 2);
				}
			}
		}
		else
		{
			//cross product with d-ntmp gammas and sigma1
			for(k=0;k<(2*i+1);k++)
			{
				Multilambda(tmp2[k], tmp[k], gens[0], pow(2,i), 2);
			}
			//cross product with unity and sigma2 (and sigma 3 if ntmp=2)
			if(ntmp!=0)
			{
				tmp3 = realloc(tmp3, sizeof(doublecomplex)*pow(2,i)*pow(2,i));
				memset(tmp3, 0, sizeof(doublecomplex)*pow(2,i)*pow(2,i));
				for(k=0;k<(int)(pow(2,i));k++)
				{
					tmp3[k*(int)(pow(2,i))+k].r = 1.0;
				}
				for(j=(2*i+1);j<(2*i+1+ntmp);j++)
				{
					if(j==(d-2))
						Multilambda(tmp2[j], tmp3, gens[1], pow(2,i), 2);
					if(j==(d-1))
						Multilambda(tmp2[j], tmp3, gens[2], pow(2,i), 2);
				}
			}
		}

		//copy temporary result from matrices tmp2 -> tmp
		for(j=0;j<(2*(i-1)+3);j++)
		{
			free(tmp[j]);
		}
		tmp = (doublecomplex**) malloc(sizeof(doublecomplex*)*(2*i+(ntmp+1)));
		for(j=0;j<(2*i+ntmp+1);j++)
		{
			tmp[j] = (doublecomplex*) malloc(sizeof(doublecomplex)*matsize*matsize);
			memcpy(tmp[j], tmp2[j], sizeof(doublecomplex)*matsize*matsize);
		}
	}

	if(N==1 && n==1)
		matsize=pow(2,N);
	else
		matsize=pow(2,N+1);


	// test for correctness of generated \gamma-matrices
	test = (doublecomplex*) malloc(sizeof(doublecomplex)*matsize*matsize);
	for(i=0;i<nummat;i++)
	{
		for(j=i;j<nummat;j++)
		{
			memset(test, 0, sizeof(doublecomplex)*matsize*matsize);
			AntiComm(test, tmp[i], tmp[j], matsize, 1.0);
			if(i!=j)
			{
				check = ScanMatrixZero(test, matsize);
				if(check!=0)
				{
					printf("------------- GammaMatrices computation -------------------------------\n");
					printf("Error in Computation! Element of test[%d*nummat+%d] is not zero!\n", i,j);
					printmat("AntiComm:", test, matsize);
				}
			}
			else
			{
				check = ScanMatrixDiag(test, matsize);
				if(check!=0)
				{
					printf("------------- GammaMatrices computation -------------------------------\n");
					printf("Error in Computation! Element of test[%d*nummat+%d]) is not zero / 2!\n", i,j);
					printmat("AntiComm:", test, matsize);
				}
			}
		}
	}

	//copy result to external matrices
	for(i=0;i<nummat;i++)
	{
		memcpy(pgamma[i], tmp[i], sizeof(doublecomplex)*matsize*matsize);
	}

	free(test);
	free(tmp3);
	for(i=0;i<nummat;i++)
	{
		free(tmp[i]);
		free(tmp2[i]);
	}
	for(i=0;i<3;i++)
	{
		free(gens[i]);
	}

	return 0;
}
