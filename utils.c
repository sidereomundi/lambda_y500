/**
 * \file utils.c
 * \brief utils module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module for useful things
 */

#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

/**
 * \brief Calculates the spectroscopic likelihood
 * \param x a double, 1st value
 * \param y a double, 2nd value
 * \return a double, minimum of x and y
 *
 * Returns the minimum of the two values
 */
double minimum(double x, double y){
  if(x < y) return x;
  return y;
}

/**
 * \brief double comparison for sorting
 * \param a a double, 1st value
 * \param b a double, 2nd value
 * \return a double, minimum of x and y
 *
 * Returns -1 is a<b, 0 if a=b and 1 if a>b
 */
int compare_double (const void * a, const void * b){
	if( *(double*)a - *(double*)b > 0) return 1;
	else if ( *(double*)a - *(double*)b < 0) return -1;
	else return 0;
}

/**
 * \brief Normal law
 * \param x a double, location where to calculate the output value
 * \param mu a double, center of the Gaussian
 * \param sigma a double, width of the Gaussian
 * \return a double, value of the Gaussian at x
 *
 * Returns the value of the Normal law at the input location
 */
double normal(double x, double mu, double sigma){
	return exp(-0.5* pow((x-mu)/sigma, 2))/(sigma*sqrt(2.*M_PI));
}

/**
 * \brief ND Normal law without correlation
 * \param x a pointer to double, location where to calculate the output value
 * \param mu a pointer double, array of centers of the Gaussian
 * \param sigma a pointer to double, array widths of the Gaussian
 * \return a double, value of the Gaussian at x
 *
 * Returns the value of the ND Normal law at the input location with a diagonal covariance matrix
 */
double normalND(double* x, double* mu, double* sigma, int N){
	int i;
	double result = 1.;
	for(i = 0; i < N; i++){
	  result *= normal(x[i], mu[i], sigma[i]) * (sigma[i]*sqrt(2.*M_PI));
	}
	return result;
}

/**
 * \brief ND Normal law
 * \param x a pointer to double, location where to calculate the output value
 * \param mu a pointer double, array of centers of the Gaussian
 * \param cov a pointer to MyMatrix, covariance matrix
 * \return a double, value of the Gaussian at x
 *
 * Returns the value of the ND Normal law at the input location
 */
double normalMultivariate(double* x, double* mu, MyMatrix* cov){

  int i;
  double* xmu;
  MyMatrix* cov1;
  double* cov1_times_xmu;
  double det;
  double arg;

  /* calculate the vector x - mu */
  xmu = malloc(cov->N*sizeof(double));
  for(i = 0; i < cov->N; i++){
    xmu[i] = x[i] - mu[i];
  }

  /* inverse the covariance matrix */
  cov1 = InvMat(cov, &det);

  /* cov-1 * (x-mu) */
  cov1_times_xmu = MatVectProd(cov1, xmu);
  /* (x-mu) * cov-1 * (x-mu) */
  arg = VectVectProd(cov->N, xmu, cov1_times_xmu);

  free(xmu);
  for(i = 0; i < cov1->N; i++) free(cov1->m[i]);
  free(cov1->m);
  free(cov1);
  free(cov1_times_xmu);

  return exp(-0.5*arg)/(sqrt(det)*pow(2.*M_PI, cov->N));

}

/**
 * \brief Lognormal law
 * \param x a double, location where to calculate the output value
 * \param mu a double, center of the Gaussian
 * \param sigma a double, width of the Gaussian
 * \return a double, value of the Gaussian at x
 *
 * Returns the value of the Lognormal law at the input location
 */
double lognormal(double x, double mu, double sigma){
	return exp(-0.5* pow((log(x)-log(mu))/sigma, 2.))/(x*sigma*sqrt(2.*M_PI));
}

