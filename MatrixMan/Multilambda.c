#include<f2c.h>
#include<stdio.h>

/*
 * computes tensor-product between two matrices; normally either \gamma-matrix or SU(N) generator and matrix generated in simulation;
 * used in 2MM,GammaMatrices, HMCYM, HMC,FreeGaussian
 */
void Multilambda(doublecomplex *out, doublecomplex *mat, doublecomplex *lambda, int size, int sizelambda)
{
	int i,j,k,l;

	for(i=0; i<size; i++)
	{
		for(k=0; k<sizelambda; k++)
		{
			for(j=0; j<size; j++)
			{
				for(l=0; l<sizelambda; l++)
				{
					out[(i*sizelambda+k)*(size*sizelambda) + j*sizelambda+l].r += mat[i*size+j].r * lambda[k*sizelambda+l].r -
							mat[i*size+j].i * lambda[k*sizelambda+l].i;
					out[(i*sizelambda+k)*(size*sizelambda) + j*sizelambda+l].i += mat[i*size+j].r * lambda[k*sizelambda+l].i +
							mat[i*size+j].i * lambda[k*sizelambda+l].r;
				}
			}
		}
	}
}
