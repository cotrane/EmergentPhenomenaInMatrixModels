#include <stdio.h>
#include <stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>


int GenStates(int **ppstate, int **ppopladder, int **ppMultiplet, int N, int dim)
{
	int i,j=1,k,l,n,m,check,st,st1,op1=0,op2=0;
	int count=0, count1=0;
	int *tempst;

	tempst = (int*) malloc(sizeof(int)*N);
	memset(tempst, 0, N*sizeof(int));

	//generating states
	i=1;
	st=1;
	st1=0;
	start2:
	start1:	if(i>((N*N-1)-(N-1)))
		{
			i=1;
			st1++;
			count1=0;
		}
		if(st1<dim)
		{
			count=0;
			for(j=1;j<N;j++)
			{
				for(k=0;k<j;k++)
				{
					if(count==count1)
					{
						op1=k;
						op2=j;
						goto exit;
					}
					count++;
				}
			}
			exit: count1++;
			if(*(ppstate[st1]+op1)>0)
			{
				tempst[op1] = *(ppstate[st1]+op1) + *(ppopladder[i]+op1);
				tempst[op2] = *(ppstate[st1]+op2) + *(ppopladder[i]+op2);
				//change 3rd value for operator
				for(n=0;n<N;n++)
				{
					if(n!=op1 && n!=op2)
					{
						tempst[n] = *(ppstate[st1]+n);
					}
				}
				//check if state already exists
				for(l=0;l<dim;l++)
				{
					if(l!=st)
					{
						check=0;
						for(m=0;m<N;m++)
						{
							if(*(ppstate[l]+m)==tempst[m])
							{
								check+=1;
							}
						}
						if(check==N)
						{
							*(ppMultiplet[i]+st1)=l;
							i+=2;
							goto start2;
						}
					}
				}
				//if not, state is copied to pst here
				for(n=0;n<N;n++)
				{
					*(ppstate[st]+n) = tempst[n];
				}
				*(ppMultiplet[i]+st1)=st;
				if(st<(dim-1))
				{
					st++;
				}
				i+=2;
				goto start1;
			}
			else
			{
				*(ppMultiplet[i]+st1)=st1;
				i+=2;
				goto start1;
			}
		}


	return 0;
}
