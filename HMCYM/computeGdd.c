/*
 * computeGdd.c
 *
 *  Created on: 22 Oct 2013
 *      Author: tkaltenbrunner
 */

#include<stdlib.h>
#include<f2c.h>
#include<clapack.h>
#include"MatrixMan.h"

/*
 * compute stress-energy tensor a la Nishimura; Gdd = Tr(X_{\mu}X_{\nu}) / N
 */
int computeGdd(doublecomplex *out, doublecomplex **in, int nummat, int size)
{
	int i,j;
	doublecomplex val;

	for(i=0;i<nummat;i++){
		for(j=i;j<nummat;j++){
			val.r=0; val.i=0;
			diagMultiCplx(&val, in[i], in[j], size, 1.0);
			if(i!=j){
				out[i*nummat+j].r = val.r;
				out[i*nummat+j].i = val.i;

				out[j*nummat+i].r = val.r;
				out[j*nummat+i].i = -val.i;
			}
			else{
				out[i*nummat+i].r = val.r;
				out[i*nummat+i].i = val.i;
			}
		}
	}

	return 0;
}
