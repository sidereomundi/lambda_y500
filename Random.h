/**
 * \file Random.h
 * \brief Random module header
 * \date on: 9 September 2011
 * \author Gurvan
 *
 * Module for random number generation
 */

#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "linalg.h"

#ifndef RANDOM_H_
#define RANDOM_H_

void InitRandom(const gsl_rng_type** T,
		gsl_rng** r);
void FreeRandom(const gsl_rng_type** T,
		gsl_rng** r);
double GetNormal(const gsl_rng_type * T,
		 gsl_rng * r, 
		 double mu, 
		 double sigma);
double GetUniform(const gsl_rng_type * T,
		  gsl_rng * r, 
		  double min, 
		  double max);
void GetNormalMultivariate(const gsl_rng_type * T,
			   gsl_rng * r,
			   int Nrand, 
			   double* mu, 
			   MyMatrix* cov, 
			   double** out,
			   double* bounds_down, 
			   double* bounds_up);


#endif /* RANDOM_H_ */
