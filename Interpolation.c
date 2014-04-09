/**
 * \file Interpolation.c
 * \brief Interpolation module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to interpolate functions
 */

#include "stdlib.h"
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "utils.h"
#include "Interpolation.h"

/**
 * \brief Interpolates a 1D function at a given location
 * \param x a pointer to double, the x axis values
 * \param vect a pointer to double, the y(x) axis values
 * \param len_x an int, the length of the vect array
 * \param x0 a double, where to interpolate
 * \return a double, the value of the interpolated value
 *
 * Does a linear interpolation of a 1D function at a given location
 */
double Interp1D_Scalar(double* x, double* y, int len_x, double x0){

    double result = 0.;

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_interp *interp = gsl_interp_alloc (gsl_interp_linear, len_x);
  gsl_interp_init (interp, x, y, len_x);

//    for(i=0;i<len_x;i++)
  //      printf("x[%d] %e y %e\n",i,x[i],y[i]);
    
    if (x0 < x[0])
        result = y[0];
    else if(x0 > x[len_x-1])
        result = y[len_x-1];
    else
        result = gsl_interp_eval (interp, x, y, x0, acc);
    

  gsl_interp_free (interp);
  gsl_interp_accel_free (acc);

  return result;
}


/**
 * \brief Interpolates a 1D function at a given location
 * \param x a pointer to double, the x axis values
 * \param vect a pointer to double, the y(x) axis values
 * \param len_x an int, the length of the vect array
 * \param x0 a double, where to interpolate
 * \return a double, the value of the interpolated value
 *
 * Does a spline-based interpolation of a 1D function at a given location
 */
double Interp1D_Scalar_Spline(double* x, double* y, int len_x, double x0){

  double result = 0.;

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, len_x);
  gsl_spline_init (spline, x, y, len_x);

  if (x0 < x[0]){
    result = y[0];
  }else if(x0 > x[len_x-1]){
    result = y[len_x-1];
  }else{
    result = gsl_spline_eval (spline, x0, acc);
  }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return result;
}

/**
 * \brief Interpolates a 1D function for a given vector of locations
 * \param x a pointer to double, the x axis values
 * \param y a pointer to double, the y(x) axis values
 * \param len_x an int, the length of the vect array
 * \param x0 a pointer to double, array of locations where to interpolate
 * \param y0 a pointer to double, array of interpolated values
 *
 * Does a linear interpolation of a 1D function for a given vector of locations
 */
void Interp1D_Vector(double* x, double* y, int len_x, double* x0, double* y0, int len_x0){
  int i,k;
    int j=0;
    
/*    // Count number of equal x
    for(i=0;i<len_x-1;i++)
        if(x[i+1]<=x[i])
            j++;
    // Redistribute them
    for(i=0;i<len_x-1;i++)
        if(x[i+1]<=x[i])
            for(k=0;k<j;k++)
                x[i+k] -= (j-k)*1.e-5;
*/
/*    FILE *fp;
    fp=fopen("monoton","w");
    for(i=0;i<len_x;i++)
        fprintf(fp,"%e %e\n",x[i],y[i]);
    fclose(fp);
    exit(0);
*/                                  
    
    
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_interp *interp = gsl_interp_alloc (gsl_interp_linear, len_x);
  gsl_interp_init (interp, x, y, len_x);

    
  //FILE* check = fopen("check.dat", "w");
  for(i = 0; i < len_x0; i++)
  {
    if ( x0[i] < x[0])
      y0[i] = y[0];
    else if(x0[i] > x[len_x-1])
      y0[i] = y[len_x-1];
    else
        y0[i] = gsl_interp_eval (interp, x, y, x0[i], acc);
    
    //fprintf(check, "%lf %e\n", x0[i], y0[i]);
  }
  //fclose(check);

    gsl_interp_free (interp);
    gsl_interp_accel_free (acc);

}

/**
 * \brief Interpolates a 1D function for a given vector of locations
 * \param x a pointer to double, the x axis values
 * \param y a pointer to double, the y(x) axis values
 * \param len_x an int, the length of the vect array
 * \param x0 a pointer to double, array of locations where to interpolate
 * \param y0 a pointer to double, array of interpolated values
 *
 * Does a spline-based interpolation of a 1D function for a given vector of locations
 */
void Interp1D_Vector_Spline(double* x, double* y, int len_x, double* x0, double* y0, int len_x0){
  int i;

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, len_x);
  gsl_spline_init (spline, x, y, len_x);

  //FILE* check = fopen("check.dat", "w");
  for(i = 0; i < len_x0; i++){
    if ( x0[i] < x[0]){
      y0[i] = y[0];
    }else if(x0[i] > x[len_x-1]){
      y0[i] = y[len_x-1];
    }else{
      y0[i] = gsl_spline_eval (spline, x0[i], acc);
    }
    //fprintf(check, "%lf %e\n", x0[i], y0[i]);
  }
  //fclose(check);

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

}

/**
 * \brief Interpolates a 2D function at a given location (x,y)
 * \param x a pointer to double, the x axis values
 * \param y a pointer to double, the y axis values
 * \param z a double pointer to double, the z(x,y) axis values
 * \param len_x an int, the length of the x array
 * \param len_y an int, the length of the y array
 * \param x0 a double, first coordinate for the interpolation
 * \param y0 a double, second coordinate for the interpolation
 * \return a double, interpolated velue
 *
 * Does a sequencial linear interpolation on both axis of a 2D function at a given location (x,y)
 */
double Interp2D_Scalar(double* x, double* y, double** z, int len_x, int len_y, double x0, double y0){

  int i;
  double interp_zy[len_x];
  double result = 0.;

  /* Interpolate on y */
  for (i = 0; i < len_x; i++){
    interp_zy[i] = Interp1D_Scalar(y, z[i], len_y, y0);
  }

  /* Interpolate on x */
  result = Interp1D_Scalar(x, interp_zy, len_x, x0);

  return result;
}

