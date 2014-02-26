#include <stdio.h>
#include <stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>

int Ops(int **ppopladder, int N)
{
	int j,k,num=0;
									//generating the operators
	for(j=1;j<N;j++)
	{
		for(k=0;k<j;k++)
		{
			*(ppopladder[num]+k)=1;
			*(ppopladder[num]+j)=-1;
			*(ppopladder[num+1]+k)=-1;
			*(ppopladder[num+1]+j)=1;
			num+=2;
		}
	}

//	printf("The operators are:\n");
//	for(i=0;i<((N*N-1)-(N-1));i++)
//	{
//		for(j=0;j<N;j++)
//		{
//			if(j<(N-1))
//			{
//				printf("%d ", *(ppopladder[i]+j));
//			}
//			else
//			{
//				printf("%d\n", *(ppopladder[i]+j));
//			}
//		}
//	}

	return 0;
}
