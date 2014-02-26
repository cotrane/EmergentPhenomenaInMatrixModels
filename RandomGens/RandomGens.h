/*
 * RandomGens.h
 *
 *  Created on: 25 Jan 2013
 *      Author: tkaltenbrunner
 */

#ifndef RANDOMGENS_H_
#define RANDOMGENS_H_

#include<f2c.h>
#include<prng.h>
#include<clapack.h>
#include<mtwist.h>

double gauss_randomnr(int n);
double gauss_randvarmean(int n, double var, double mean);
doublecomplex gauss_randcomplex();


#endif /* RANDOMGENS_H_ */
