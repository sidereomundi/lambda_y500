/**
 * \file Random.c
 * \brief Random module
 * \date on: 9 September 2011
 * \author Gurvan
 *
 * Module for random number generation
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>

#include "Random.h"

/**
 * \brief Inits the random number generator
 * \param T a double pointer to gsl_rng_type, gsl parameter
 * \param r a double pointer to gsl_rng, gsl parameter
 *
 * Inits the random number generator using the GSL
 */
void InitRandom(const gsl_rng_type ** T,
		gsl_rng ** r){
	/* create a generator chosen by the
	environment variable GSL_RNG_TYPE */
	long seed;

	gsl_rng_env_setup();
	*T = gsl_rng_default;
	*r = gsl_rng_alloc (*T);

	seed = time (NULL) * getpid();
	gsl_rng_set (*r, seed);

}

/**
 * \brief Deallocates the random number generator
 * \param T a double pointer to gsl_rng_type, gsl parameter
 * \param r a double pointer to gsl_rng, gsl parameter
 *
 * Deallocates the GSL random number generator
 */
void FreeRandom(const gsl_rng_type ** T,
		gsl_rng ** r){
    gsl_rng_free (*r);
}

/**
 * \brief Get a random number in a Normal distribution
 * \param T a pointer to gsl_rng_type, gsl parameter
 * \param r a pointer to gsl_rng, gsl parameter
 * \param mu a double, Gaussian center
 * \param sigma a double, Gaussian width
 *
 * Get a random number in a Normal distribution given mu and sigma
 */
double GetNormal(const gsl_rng_type * T,
		 gsl_rng * r,
		 double mu, 
		 double sigma){
	double k;

	/* get a random number chosen from
	the normal distribution with
	parameters mu and sigma */
    k = mu + gsl_ran_gaussian (r, sigma);

    return k;
}

/**
 * \brief Get a random number in a Uniform distribution
 * \param T a pointer to gsl_rng_type, gsl parameter
 * \param r a pointer to gsl_rng, gsl parameter
 * \param min a double, minimum boundary
 * \param max a double, maximum boundary
 *
 * Get a random number in a Uniform distribution given boundaries
 */
double GetUniform(const gsl_rng_type * T,
		  gsl_rng * r,
		   double min,
		   double max){
	double k;

	/* get a random number chosen from
	the uniform distribution with
	parameters min and max */
    k = gsl_rng_uniform (r)*(max-min) + min;

    return k;

}

/**
 * \brief Get a random number in a Normal Multivariate distribution
 * \param T a pointer to gsl_rng_type, gsl parameter
 * \param r a pointer to gsl_rng, gsl parameter
 * \param Nrand an int, number of random points
 * \param mu a pointer to double, array of Gaussian centers
 * \param cov a pointer to MyMatrix, covariance matrix
 * \param out a double pointer to double, array of generated coordinates
 * \param bounds_down a pointer to double, lower boundary array
 * \param bounds_up a pointer to double, upper boundary array
 *
 * Get a random number in a Normal Multivariate distribution given a center vector and a covariance matrix.
 */
void GetNormalMultivariate(const gsl_rng_type * T,
			   gsl_rng * r,
			   int Nrand, 
			   double* mu, 
			   MyMatrix* cov, 
			   double** out,
			   double* bounds_down, 
			   double* bounds_up){

  int i, j;
  MyMatrix* eigvec;
  double* eigval;
  double* randvec;
  double* randvecrot;
  short outofrange;

  /* diagonalize the covariance matrix */
  eigvec = malloc(sizeof(MyMatrix));
  eigvec->N = cov->N;
  eigvec->m = malloc(eigvec->N*sizeof(double*));
  for(i = 0; i < cov->N; i++) eigvec->m[i] = malloc(eigvec->N*sizeof(double));
  eigval = malloc(cov->N*sizeof(double));

  Diagonalize(cov, eigvec, eigval);
  
  /*
  printf("eigenvalues\n");
  for(j = 0; j < cov->N; j++) printf("%e\t", eigval[j]);
  printf("\n");
  printf("eigenvectors\n");
  for(i = 0; i < cov->N; i++){
    for(j = 0; j < cov->N; j++)
      printf("%lf\t", eigvec->m[i][j]);
    printf("\n");
  }  
  */

  /* Generate normal distributions along the eigen vectors */
  randvec = malloc(cov->N*sizeof(double));
  i = 0;
  while(i<Nrand){
    for(j = 0; j < cov->N; j++){
      if(eigval[j]<0.){
	printf("negative %d %lf\n", j, eigval[j]);
	exit(0);
      }else if(isnan(eigval[j])){
 	printf("nan %d %lf\n", j, eigval[j]);
	exit(0);
     }else if(isinf(eigval[j])){
	printf("inf %d %lf\n", j, eigval[j]);
	exit(0);
      }
      randvec[j] = GetNormal(T, r, 0., sqrt(eigval[j]));
    }
    randvecrot = MatVectProd(eigvec, randvec);

    /* test if the random vector is in boundaries */
    outofrange = 0;
    for(j = 0; j < cov->N; j++){
      if(randvecrot[j] + mu[j] < bounds_down[j] || randvecrot[j] + mu[j] > bounds_up[j]){
	outofrange = 1;
      }
    }

    /* if so, fill in the output */
    if (!outofrange){
      for(j = 0; j < cov->N; j++){
	out[i][j] = randvecrot[j] + mu[j];
      }
      i++;
    } /*else{
      printf("out of range\n");
      for(j = 0; j < cov->N; j++){
      	printf("%lf\t%lf\n", randvec[j], randvecrot[j] + mu[j]);
      }
      }*/
    
    free(randvecrot);
  }

  /* free */
  for(i = 0; i < eigvec->N; i++) free(eigvec->m[i]);
  free(eigvec->m);
  free(eigvec);
  free(eigval);
  free(randvec);

}


