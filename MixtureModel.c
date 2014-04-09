/**
 * \file MixtureModel.c
 * \brief MixtureModel module
 * \date on: May 2011
 * \author Gurvan
 *
 * Module to handle a mixture model for PMC
 */

#include <stdio.h>
#include <stdlib.h>

#include "MixtureModel.h"
#include "Random.h"
#include "utils.h"


/**
 * \brief Inits mixture model
 * \param T a pointer to gsl_rng_type, for random number generation
 * \param r a pointer to gsl_rng, for random number generation
 * \param Ncomp an int, number of components
 * \param Nparam an int, number of parameters
 * \param priorcenter a pointer to double, center of input Gaussian priors
 * \param priorwidth a pointer to double, width of input Gaussian priors
 * \return a pointer to mix_model, initialized mixture model
 *
 * Allocates and inits mixture model using random values from priors given the number of components and number of parameters
 */
mix_model* InitMixModel(const gsl_rng_type* T,
			gsl_rng* r,
			int Ncomp, 
			int Nparam, 
			double* priorcenter,
			double* priorwidth){

  int i, j, k;
  mix_model* mixmodel;
  
  /* memory allocation */
  mixmodel = malloc(sizeof(mix_model));
  mixmodel->Ncomp = Ncomp;
  mixmodel->cov_vect = malloc(Ncomp*sizeof(MyMatrix*));
  mixmodel->mu_vect = malloc(Ncomp*sizeof(double*));
  mixmodel->alpha_vect = malloc(Ncomp*sizeof(double));
  mixmodel->validcomp = malloc(Ncomp*sizeof(short));
  for(i = 0; i < Ncomp; i++){
    mixmodel->cov_vect[i] = malloc(sizeof(MyMatrix));
    mixmodel->cov_vect[i]->N = Nparam;
    mixmodel->cov_vect[i]->m = malloc(Nparam*sizeof(double*));
    for(j = 0; j < Nparam; j++){
      mixmodel->cov_vect[i]->m[j] = malloc(Nparam*sizeof(double));
    }
    mixmodel->mu_vect[i] = malloc(Nparam*sizeof(double));
  }

  /* initialize the mixture model parameters */
  for(i = 0; i < Ncomp; i++){
    mixmodel->validcomp[i] = 1;
    mixmodel->alpha_vect[i] = 1. / (double) (Ncomp);
    for(j = 0; j < Nparam; j++){
      //mixmodel->mu_vect[i][j] = 1. + GetNormal(0., 2.);
      mixmodel->mu_vect[i][j] = GetNormal(T, r, priorcenter[j], priorwidth[j]);
      for(k = j; k < Nparam; k++){
	//mixmodel->cov_vect[i]->m[j][k] = ( j == k ? GetNormal(5, 1.) : GetNormal(1., 1.) );
	mixmodel->cov_vect[i]->m[j][k] = ( j == k ? GetNormal(T, r, priorwidth[j]*priorwidth[j], priorwidth[j]*priorwidth[j]*0.1) : GetNormal(T, r, 0., priorwidth[j]*priorwidth[k]*1.e-4) );
      }
      for(k = 0; k < j; k++){
	mixmodel->cov_vect[i]->m[j][k] = mixmodel->cov_vect[i]->m[k][j];
      }
    }
  }

  for(i = 0; i < Ncomp; i++){
    printf("mixmodel %d\n", i);
    for(j = 0; j < Nparam; j++){
      printf("%lf\t", mixmodel->mu_vect[i][j] );
    }
    printf("\n");
    for(j = 0; j < Nparam; j++){
      for(k = j; k < Nparam; k++){
	printf("%lf\t", mixmodel->cov_vect[i]->m[j][k]);
      }
      printf("\n");
    }
  }
  printf("\n");
 

  return mixmodel;

}

/**
 * \brief Inits mixture model from an input file
 * \param filename a pointer to char, input file path
 * \return a pointer to mix_model, initialized mixture model
 *
 * Allocates and inits mixture model from an input file, for instance from a previous PMC chain
 */
mix_model* SetMixModel(char* filename){

  int Ncomp, Nparam;
  int i, j, k, n;
  mix_model* mixmodel;
  char buffer[8192];
  char* b;
  int offset;
  FILE* myfile;
  double tmp;

  /* read the input file */
  myfile = fopen(filename, "r");
  fgets(buffer, sizeof(buffer), myfile);
  n = sscanf(buffer, "%d %d", &Ncomp, &Nparam);

  /* memory allocation */
  mixmodel = malloc(sizeof(mix_model));
  mixmodel->Ncomp = Ncomp;
  mixmodel->cov_vect = malloc(Ncomp*sizeof(MyMatrix*));
  mixmodel->mu_vect = malloc(Ncomp*sizeof(double*));
  mixmodel->alpha_vect = malloc(Ncomp*sizeof(double));
  mixmodel->validcomp = malloc(Ncomp*sizeof(short));
  for(i = 0; i < Ncomp; i++){
    mixmodel->cov_vect[i] = malloc(sizeof(MyMatrix));
    mixmodel->cov_vect[i]->N = Nparam;
    mixmodel->cov_vect[i]->m = malloc(Nparam*sizeof(double*));
    for(j = 0; j < Nparam; j++){
      mixmodel->cov_vect[i]->m[j] = malloc(Nparam*sizeof(double));
    }
    mixmodel->mu_vect[i] = malloc(Nparam*sizeof(double));
  }

  /* initialize the mixture model parameters */
  fgets(buffer, sizeof(buffer), myfile);
  b = buffer;
  for(i = 0; i < Ncomp; i++){
    sscanf(b, "%lf %n", &tmp, &offset);
    b += offset;
    mixmodel->validcomp[i] = 1;
    mixmodel->alpha_vect[i] = tmp;
  }
  for(i = 0; i < Ncomp; i++){
    fgets(buffer, sizeof(buffer), myfile);
    b = buffer;
    for(j = 0; j < Nparam; j++){
      sscanf(b, "%lf %n", &tmp, &offset);
      b += offset;
      mixmodel->mu_vect[i][j] = tmp;
    }
    for(j = 0; j < Nparam; j++){
      fgets(buffer, sizeof(buffer), myfile);
      b = buffer;
      for(k = j; k < Nparam; k++){
	sscanf(b, "%lf %n", &tmp, &offset);
	b += offset;
	mixmodel->cov_vect[i]->m[j][k] = tmp;
      }
      for(k = 0; k < j; k++){
	mixmodel->cov_vect[i]->m[j][k] = mixmodel->cov_vect[i]->m[k][j];
      }
    }
  }
  fclose(myfile);

  for(i = 0; i < Ncomp; i++){
    printf("mixmodel %d\t%lf\n", i, mixmodel->alpha_vect[i]);
    for(j = 0; j < Nparam; j++){
      printf("%lf\t", mixmodel->mu_vect[i][j] );
    }
    printf("\n");
    for(j = 0; j < Nparam; j++){
      for(k = j; k < Nparam; k++){
	printf("%lf\t", mixmodel->cov_vect[i]->m[j][k]);
      }
      printf("\n");
    }
  }
  printf("\n");
 
  return mixmodel;

}

/**
 * \brief Deallocates a mixture model
 * \param a pointer to mix_model, mixture model structure
 *
 * Deallocates a mixture model structure
 */
void FreeMixModel(mix_model* mixmodel){

  int i, j;

  for(i = 0; i < mixmodel->Ncomp; i++){
    for(j = 0; j < mixmodel->cov_vect[i]->N; j++){
      free(mixmodel->cov_vect[i]->m[j]);
    }
    free(mixmodel->cov_vect[i]->m);
    free(mixmodel->cov_vect[i]);
    free(mixmodel->mu_vect[i]);
  }
  free(mixmodel->alpha_vect);
  free(mixmodel->validcomp);
  free(mixmodel);

}

/**
 * \brief Updates a mixture model
 * \param a pointer to mix_model, mixture model structure
 * \param N an int, number of random realizations in the PMC iteration
 * \param x a double pointer to double, array or random values for parameters
 * \param pi a pointer to double, output array of importance (importance sampling)
 * \param perplexity a pointer to double, output perplexity of the PMC (covergence criterion)
 * \param ESS a pointer to double, output Effective Sample Size (covergence criterion)
 *
 * Updates a mixture model structure given a PMC iteration
 */
void UpdateMixModel(mix_model* mixmodel, int N, double** x, double* pi, double* perplexity, double* ESS){
  
  int i, j, k, l;
  double** rho;
  double* pp;
  double* wx;
  double sumwx;
  double H;
  int n;

  /* calculate the mixture model at each point */
  rho = malloc(mixmodel->Ncomp*sizeof(double*));
  for(i = 0; i < mixmodel->Ncomp; i++) rho[i] = malloc(N*sizeof(double));
  pp = malloc(N*sizeof(double));
  wx = malloc(N*sizeof(double));

  /* calculate rho, weight, posterior pi and mixture model */
  sumwx = 0.;
  for(i = 0; i < N; i++){
    pp[i] = 0.;
    for(j = 0; j < mixmodel->Ncomp; j++){
      if(mixmodel->validcomp[j]){
	rho[j][i] = mixmodel->alpha_vect[j] * normalMultivariate(x[i], mixmodel->mu_vect[j], mixmodel->cov_vect[j]);
	pp[i] += rho[j][i];
      }else{
	rho[j][i] = 0.;
      }
    }
    for(j = 0; j < mixmodel->Ncomp; j++){
      if(mixmodel->validcomp[j]){
	rho[j][i] /= pp[i];
      }
    }
    wx[i] = pi[i] / pp[i];
    sumwx += wx[i];

  }
  for(i = 0; i < N; i++){
    wx[i] /= sumwx;
  }

  /* update parameters */
  for(j = 0; j < mixmodel->Ncomp; j++){
    mixmodel->alpha_vect[j] = 0.;
    if(mixmodel->validcomp[j]){
      for(i = 0; i < N; i++){
	mixmodel->alpha_vect[j] += wx[i] * rho[j][i];
      }
      if(mixmodel->alpha_vect[j] < 0.01){
	mixmodel->validcomp[j] = 0;
	printf("Mixture model component %d discarded\n", j);
      }
    }
  }

  for(j = 0; j < mixmodel->Ncomp; j++){
    if (!mixmodel->validcomp[j]) continue;
    for(k = 0; k < mixmodel->cov_vect[j]->N; k++){
      mixmodel->mu_vect[j][k] = 0.;
      for(i = 0; i < N; i++){
	mixmodel->mu_vect[j][k] += wx[i] * x[i][k] *rho[j][i] / mixmodel->alpha_vect[j];
      }
      
      for(l = 0; l < mixmodel->cov_vect[j]->N; l++){
	mixmodel->cov_vect[j]->m[k][l] = 0.;
	for(i = 0; i < N; i++){
	  mixmodel->cov_vect[j]->m[k][l] +=  wx[i] * (x[i][k] - mixmodel->mu_vect[j][k]) * (x[i][l] - mixmodel->mu_vect[j][l]) * rho[j][i] / mixmodel->alpha_vect[j];
	}//i
      }//l	    

    }//k
  }//j

  for(i = 0; i < mixmodel->Ncomp; i++){
    if (!mixmodel->validcomp[i]) continue;
    printf("mixmodel update %d\n", i);
    for(j = 0; j < mixmodel->cov_vect[i]->N; j++){
      printf("%lf\t", mixmodel->mu_vect[i][j] );
    }
    printf("\n");
    for(j = 0; j < mixmodel->cov_vect[i]->N; j++){
      for(k = j; k < mixmodel->cov_vect[i]->N; k++){
	printf("%lf\t", mixmodel->cov_vect[i]->m[j][k]);
      }
      printf("\n");
    }
  }
  printf("\n");

  /* Shannon entropy */
  H = 0.;
  n = 0;
  for(i = 0; i < N; i++){
    if(wx[i]>0.){
      H -= wx[i] * log(wx[i]);
      n++;
    }
  }
  /* Normalized perplexity */
  *perplexity = exp(H) / (double) n;

  /* Effective sample size */
  *ESS = 0.;
  for(i = 0; i < N; i++){
    if(wx[i]>0.){
      *ESS += wx[i]*wx[i];
    }
  }
  *ESS = 1. / *ESS;

  for(i = 0; i < mixmodel->Ncomp; i++) free(rho[i]);
  free(rho);
  free(pp);
  free(wx);

}



