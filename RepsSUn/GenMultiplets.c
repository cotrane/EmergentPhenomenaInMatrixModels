#include <stdio.h>
#include <stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>

//generating the Multiplets

int GenMultiplets(int **ppMultiplet, int **ppMultinum, int N, int dim, int numMultis)
{
	int i,j,k,l,m;

	//counts operators
	for(i=1;i<((N*N-1)-(N-1));i+=2)
	{
		//numbers the Multiplet for one operator
		l=0;
		//counts states
		for(j=0;j<dim;j++)
		{
			//checks if state appears in different Multiplet; we just want to start from highest weight in this multiplet!
			for(k=0;k<dim;k++)
			{
				if(k!=j && *(ppMultiplet[i]+k)==j)
				{
					goto start3;
				}
			}
			m=j;
			if(*(ppMultinum[i]+l)==0)
			{
				*(ppMultinum[i]+l)=1;
			}
			//runs through Multiplet and counts states
			while(*(ppMultiplet[i]+m)!=m)
			{
				*(ppMultinum[i]+l)+=1;
				m=*(ppMultiplet[i]+m);
			}
			l++;
			start3:;
		}
	}


	return 0;
}

