/**
 * \file Statistics.c
 * \brief Statistics module
 * \date on: 11 August 2011
 * \author Gurvan
 *
 * Module to calculate the likelihood
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Statistics.h"
#include "Integration.h"
#include "Interpolation.h"
#include "utils.h"


/**
 * \brief Inits the tabulated error function
 * \param filename a pointer to char, path of the input file
 * \return a pointer to erf_struct, output error function structure
 *
 * Allocates and initializes the tabulated error function
 */
erf_struct* InitErf(char* filename){

  FILE* myfile;
  int i;
  long double x, y;
  int N = 27999;
  erf_struct* errfunc;
 
  /* allocate memory */
  errfunc = malloc(sizeof(erf_struct));
  errfunc->x = malloc(N * sizeof(long double));
  errfunc->y = malloc(N * sizeof(long double));
  
  /* read the data file */
  myfile = fopen(filename, "r");
  i = 0;
  for (i = 0; i < N; i++){
    fscanf(myfile, "%Lf", &x);
    fscanf(myfile, "%Lf", &y);
    errfunc->x[i] = x;
    errfunc->y[i] = y*2.;
  }
  fclose(myfile);

  return errfunc;
}

/**
 * \brief Deallocates the tabulated error function
 * \param a double pointer to erf_struct, input error function structure
 *
 * Deallocates the tabulated error function
 */
void FreeErf(erf_struct** errfunc){
  free((*errfunc)->x);
  free((*errfunc)->y);
  free(*errfunc);
}

/**
 * \brief Inverse error function
 * \param errfunc a pointer to erf_struct, tabulated error function
 * \param y0 a long double, the input value
 *
 * Calculates the inverse error function from the tabulated error function
 */
double inverf(erf_struct* errfunc, long double y0){

  int N = 27999;
  int i;
  int ibreak = 0;

  /* look for the location of y0 */
  if (y0 < errfunc->y[0]) return errfunc->x[0];
  if (y0 > errfunc->y[N-2]) return errfunc->x[N-1];
  
  for(i = 0; i < N; i++){
    if(errfunc->y[i] > y0){
      ibreak = i-1;
      break;
    }
  }
  //printf("ibreak = %d\t%Lf\t%Lf\t%Lf\n", ibreak, errfunc->y[ibreak], errfunc->y[i], y0);
  return (double) (errfunc->x[ibreak] + (y0 - errfunc->y[ibreak]) * (errfunc->x[ibreak+1]-errfunc->x[ibreak])/(errfunc->y[ibreak+1]-errfunc->y[ibreak]));
}


/**
 * \brief Get the pull in linear scale
 * \param f a pointer to Function_1D, input distribution 
 * \param meas_value a double, input measured value
 * \return a double, value of the pull
 *
 * Calculates the pull making the approximation of a symetrc distribution in linear scale
 */
double Get_Pull(Function_1D* f, double meas_value){
	int i;
	double left;
	double right;
	double middle;
	double sigma;
	double frac = 0.16;

	/* calculate the location of the left 1sigma value */
	i = GetLocation_IntTrap_Left(f->x, f->y, f->N, frac);
	if(i >= f->N || i == 0){
		printf( "Error: Statistics::Get_Pull: left index out of bound, return 0\n");
		return 0;
	}
	left = f->x[i];
	/* calculate the location of the right 1sigma value */
	i = GetLocation_IntTrap_Right(f->x, f->y, f->N, frac);
	if(i >= f->N || i == 0){
		printf( "Error: Statistics::Get_Pull: right index out of bound, return 0\n");
		return 0;
	}
	right = f->x[i];
	middle = (right + left) / 2.;
	sigma = (right - left) / 2.;
	//printf("%lf,\t%lf\t%lf\n", left, middle, right);

	return (meas_value - middle) / sigma;
}

/**
 * \brief Get the pull
 * \param f a pointer to Function_1D, input distribution
 * \param meas_value a double, input measured value
 * \param errfunc a pointer to erf_struct, tabulated error function
 * \return a double, value of the pull
 *
 * Calculates the pull without approximation
 */
double Get_GaussPull(Function_1D* f, double meas_value, erf_struct* errfunc){
	int i, n;
	double dx, sum;
	double* y;
	double result;

	/* find where is the measured value */
	i = 0;
	while(1){
		if(f->x[i] > meas_value || i == f->N-1) break;
		i++;
	}

	/* if out of bound */
	if(i == f->N-1 || i == 0){
	  printf( "Error: Statistics::Get_GaussPull: value %lf out of bounds [%lf, %lf]\n", meas_value, f->x[0], f->x[f->N-1]);
		return WRONG_PULL;
	}

	/* Integrate from -inf to meas_value */
	n = i;
	if(n-1 > 0)
    {
        y = (double*) malloc(n * sizeof(double));
        for(i=0; i < n; i++)
            y[i] = f->y[i];
        
        dx = f->x[1] - f->x[0];
        sum = IntTrap_conststep(dx, y, n);
        free(y);
        
        if (isnan(sum)){
            printf("NANSUM before addition\n");
            return -10000.;
        }        
	}
    else
        sum = 0.;
	n--;
    // Interpolate to meas_value
    sum += (meas_value - f->x[n]) * (f->y[n] + (meas_value - f->x[n]) * (f->y[n+1] - f->y[n]) / (f->x[n+1] - f->x[n]));
	
    if (isnan(sum)) printf("NANSUM %d\t%lf\t%lf\t%lf\t%lf\t%lf\n", n, f->x[n], f->x[n+1], meas_value, f->y[n], f->y[n+1]);
	/* return the corresponding number of sigma for a gaussian */
	/* cdf = 0.5 * (1+erf((x-mu)/sqrt(2))) */
	result = inverf(errfunc, 2.*sum - 1.);
	//if(isnan(sum)) printf("nan result: %lf %d %lf %lf %lf %lf %lf %lf\n", sum, n, dx, meas_value, f->x[n], f->x[n+1], f->y[n], f->y[n+1]);
	//printf("gausspull %lf\t%lf\n", result, sum);
	if(isnan(result)) return WRONG_PULL;
	return result;
}

/**
 * \brief Get the pull in log scale
 * \param f a pointer to Function_1D, input distribution 
 * \param meas_value a double, input measured value
 * \return a double, value of the pull
 *
 * Calculates the pull making the approximation of a symetrc distribution in log scale
 */
double Get_LogPull(Function_1D* f, double meas_value){
	int i;
	double left;
	double right;
	double middle;
	double sigma;
	double frac = 0.16;
	double logx[f->N];
	double logx_conststep[f->N];
	double step;
	double interpy_logx_conststep[f->N];
	double norm;

	/* Interpolate over a constant log step grid */
	step = (log(f->x[f->N-1]) - log(f->x[0]))/ ((double) f->N);
	for(i = 0; i < f->N; i++){
		logx[i] = log(f->x[i]);
		logx_conststep[i] = log(f->x[0])+i*step;
	}
	Interp1D_Vector(logx, f->y, f->N, logx_conststep, interpy_logx_conststep, f->N);
	/* normalize the distribution is log space */
	norm = IntTrap_conststep(step, interpy_logx_conststep, f->N);
	for(i = 0; i < f->N; i++) interpy_logx_conststep[i] /= norm;

	/* calculate the location of the left 1sigma value */
	i = GetLocation_IntTrap_Left(logx_conststep, interpy_logx_conststep, f->N, frac);
	if(i >= f->N || i == 0){
		printf("Error: Statistics::Get_LogPull: left index out of bound\n");
		return WRONG_PULL;
	}
	left = logx_conststep[i];
	/* calculate the location of the right 1sigma value */
	i = GetLocation_IntTrap_Right(logx_conststep, interpy_logx_conststep, f->N, frac);
	if(i >= f->N || i == 0){
		printf( "Error: Statistics::Get_LogPull: right index out of bound\n");
		return WRONG_PULL;
	}
	right = logx_conststep[i];
	middle = (right + left) / 2.;
	sigma = (right - left) / 2.;
	//printf("%lf,\t%lf\t%lf\n", left, middle, right);

	return (log(meas_value) - middle) / sigma;
}

/**
 * \brief KS test intermediate function
 * \param alam a double
 * \return a double
 *
 * from Numerical Recipes
 */
double Get_KS_p(double alam){
	/* Following Numerical Recipes */
	int j;
	double a2, fac=2.0, sum=0., term, termbf=0.;

	a2 = -2. * alam*alam;
	for(j = 1; j <= 100; j++){
		term = fac*exp(a2*j*j);
		sum += term;
		//printf("%d %e %e\n", j, term, sum);
		if(fabs(term) <= EPS1*termbf || fabs(term) <= EPS2*sum) return sum;
		fac = -fac;
		termbf = fabs(term);
	}
	return 1.;
}

/**
 * \brief KS test
 * \param distrib_in a pointer to double, input distribution 
 * \param len an int, length of the array distrib_in
 * \return a double, p-value of the KS test
 *
 * Calculates the KS test p-value
 */
double KS_Pvalue(double* distrib_in, int len){
	/* Following Numerical Recipes */
	int i;
	double x;
	double cumulative;
	double normal_cumulative;
	double D;
	double absdiff;
	int N = 0;
	double* distrib;
	double result = 0.;

	/* count the number of valid values */
	/*for(i = 0; i < len; i++){
			if(distrib_in[i]>-1000) N++;
	}*/
	N = len;
	/* copy the valid values */
	distrib = (double*) malloc(N * sizeof(double));
	N = 0;
	for(i = 0; i < len; i++){
			if(distrib_in[i]>-1000.) distrib[N++] = distrib_in[i];
			else distrib[N++] = 5.; // arbitrary
	}
	len = N;

	if(len != 0){
		/* sort the distribution */
		qsort (distrib, len, sizeof(double), compare_double);

		/* compute the cumulative distributions
		 * and the maximum difference
		 */
		D = 0.;
		for(i = 0; i < len; i++){
			x = distrib[i];
			cumulative = (i+1) / ((double) len);
			normal_cumulative = 0.5*(1.+erf(x/M_SQRT2));
			absdiff = fabs(cumulative - normal_cumulative);
			if(D < fabs(absdiff)) D = absdiff;
		}

		/* return the p-value */
		len = sqrt(len);

		result = Get_KS_p( (len+0.12+0.11/len) * D );

	}

	free(distrib);

	return result;
}

