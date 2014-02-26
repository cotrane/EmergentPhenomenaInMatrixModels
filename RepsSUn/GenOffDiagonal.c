#include <stdio.h>
#include <stdlib.h>
#include<f2c.h>
#include<string.h>
#include<math.h>

/*
 * generates non-diagonal matrices of lie algebra representation
 */

int GenOffDiagonal(doublecomplex **ppmat, int **ppMultiplet, int **ppMultinum, int N, int dim, double a)
{
	int i,j,k,l,m,n,row,col,multpos,num=1,num2;
	double hw,hw1,entry;
	doublecomplex *plad;


	//generating the off-diagonal basis matrices
	for(i=1;i<((N*N-1)-(N-1));i+=2)						//counts operators
	{
		l=0;											//numbers Multiplet for one operator
		for(j=0;j<dim;j++)								//counts states
		{
			for(n=2;n<=N;n++)							//leaving the diagonal places out
			{
				num2=n*n-1;
				if(num2>num)
					break;
				else if(num==num2){
					num++;
					break;
				}
			}
			if(num>=(N*N-1))
				break;
			for(k=0;k<dim;k++)							//checks if state appears in Multiplet treated earlier
			{
				if(k!=j && *(ppMultiplet[i]+k)==j)
				{
					goto start4;
				}
			}
			m=j;
			multpos=0;
			while(*(ppMultiplet[i]+m)!=m)							//runs through Multiplet and counts states
			{
				hw = *(ppMultinum[i]+l);
				hw = (hw-1.0)/2.0;
				hw1= hw-multpos;
				entry = a/2.0 *(sqrt((hw+hw1)*(hw-hw1+1.0)));
				row=*(ppMultiplet[i]+m);
				col=m;
				if(row<col)
				{
					plad = ppmat[num-1];
					plad[row*dim+col].r=entry;
					plad[col*dim+row].r=entry;
					plad = ppmat[num];
					plad[row*dim+col].i=-entry;
					plad[col*dim+row].i=entry;
				}
				else
				{
					plad = ppmat[num-1];
					plad[row*dim+col].r=entry;
					plad[col*dim+row].r=entry;
					plad = ppmat[num];
					plad[row*dim+col].i=entry;
					plad[col*dim+row].i=-entry;
				}
				m=*(ppMultiplet[i]+m);
				multpos+=1;
			}
			l++;
			start4:;
		}
		num=num+2;
	}


	return 0;
}
