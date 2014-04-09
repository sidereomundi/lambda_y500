/**
 * \file Likelihood.h
 * \brief Likelihood module header
 * \date on: 8 September 2011
 * \author Gurvan
 *
 * Module to calculate the likelihood
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Data.h"
#include "Function.h"
#include "Statistics.h"
#include "ScalingRelation.h"
#include "Parameters.h"

#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

#ifndef DIST_H
#define DIST_H 2997.92458   // c/H*h in Mpc
#endif

#ifndef PI
#define PI 3.14159265359
#endif



double GetLikelihood_Y500cyl(Parameters_struct* params,
			     Function_1D** vect_p_m_lambda_debiased, 
			     int Ncluster, 
			     Function_1D** vect_p_Y500obs_lambda,
			     Data_structure* data );


double Likelihood_Lambdacalib(Parameters_struct* params, 
			      int Ncluster, 
			      Data_structure* data,
			      Function_2D* mass_func,
			      y0_to_Y500cyl_structure* y0_to_Y500cyl);

#endif
