/**
 * \file utils.h
 * \brief utils module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module for useful things
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>

#include "linalg.h"

/**
 * \def MIN
 * \brief minimum of two values
 */
#define MIN(x, y) x > y ? y : x

double minimum(double, double);
int compare_double (const void *, const void *);
double lognormal(double, double, double);
double normal(double, double, double);
double normalND(double*, double*, double*, int);
double normalMultivariate(double*, double*, MyMatrix*);

#endif /* UTILS_H_ */
