/**
 * \file MixtureModel.h
 * \brief MixtureModel module header
 * \date on: May 2011
 * \author Gurvan
 *
 * Module to handle a mixture model for PMC
 */

#include "linalg.h"
#include "Random.h"

#ifndef MIXMODEL_H_
#define MIXMODEL_H_

/**
 * \struct mix_model
 * \brief Mixture Model structure
 */
typedef struct{
  int Ncomp; /**< an int, Number of components */
  MyMatrix** cov_vect; /**< a double pointer to MyMatrix, array of Gaussian covariance matrices */
  double** mu_vect; /**< a double pointer to double, array or mu vectors */
  double* alpha_vect; /**< a pointer to double, array of amplitudes */
  short* validcomp; /**< a pointer to short, array of values 0 or 1. 1 for enabled component, 0 for disabled component in the mixture model */
} mix_model;

mix_model* InitMixModel(const gsl_rng_type* T,
			gsl_rng* r,
			int Ncomp, 
			int Nparam, 
			double* priorcenter,
			double* priorwidth);
void FreeMixModel(mix_model*);
mix_model* SetMixModel(char*);
void UpdateMixModel(mix_model*, int, double**, double*, double*, double*);


#endif

