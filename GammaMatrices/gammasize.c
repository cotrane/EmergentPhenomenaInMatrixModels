/*
 * gammasize.c
 *
 *  Created on: 8 Nov 2012
 *      Author: tkaltenbrunner
 */

#include<math.h>

/*
 * returns the size of \gamma-matrices corresponding to the number of matrices defined by nummat
 */
int computeGammasize(int nummat)
{
	int i,j;
	int gammasize=0;

	if(nummat>2)
	{
		//define number/size of gamma matrices; 2*N+n - dimensional
		for(i=1;i<20;i++)
		{
			for(j=1;j<4;j++)
			{
				if((nummat/(2.*i+j)) == 1.0)
				{
					if(i==1 && j==1)
					{
						gammasize=pow(2,i);
						break;
					}
					else
					{
						gammasize=pow(2,i+1);
						break;
					}
				}
			}
			if(gammasize!=0)
				break;
		}
	}
	else
		gammasize=2;

	return gammasize;
}

