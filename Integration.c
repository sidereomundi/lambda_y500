/**
 * \file Integration.c
 * \brief Integration module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to integrate a function
 */

#include <stdio.h>

#include "utils.h"
#include "Integration.h"

/**
 * \brief Integrates a function with a constant binning
 * \param h a double, size of the bin
 * \param vect a pointer to double, the function values
 * \param len_vect an int, the length of the vect array
 * \return a double, the value of the integral
 *
 * Uses the trapezoid method to integrate a function with a constant binning
 */
double IntTrap_conststep(double h, double* vect, int len_vect){
	double integral;
	int i;

	integral = 0.;
	for (i = 1; i< len_vect-1; i++) integral += vect[i];
	return h*( integral + (vect[0]+vect[len_vect-1])/2. );
}

/**
 * \brief Integrates a function with a variable binning
 * \param x_vect a pointer to double, the x axis values
 * \param y_vect a pointer to double, the y(x) axis values
 * \param len_x an int, the length of the x_vect array
 * \param len_y an int, the length of the y_vect array
 * \return a double, the value of the integral
 *
 * Uses the trapezoid method to integrate a function with a variable binning
 */
double IntTrap_varstep(double* x_vect, double* y_vect, int len_x, int len_y)
{
  double integral;
  int i;
  
  if(len_x != len_y){
    printf("Integration::IntTrap_varstep: Error, size(x) != size(y), return 0");
    return 0;
  }
  
  integral = 0.;
  for (i = 0; i< len_x-1; i++) 
    {
      integral += (y_vect[i+1]+y_vect[i]) * (x_vect[i+1]-x_vect[i]);
      //      printf("%d: intg:%f\n ",i,integral);
    }
  return integral/2.;
}

/**
 * \brief Gets a percentil location integrating from the left
 * \param x a pointer to double, the x axis values
 * \param vect a pointer to double, the y(x) axis values
 * \param len_y an int, the length of the vect array
 * \param frac a double, the given percentil
 * \return a double, the value of the percentil location
 *
 * Gets a percentil location integrating from the left
 */
int GetLocation_IntTrap_Left(double* x, double* vect, int len_vect, double frac){
	double integral;
	int i;

	integral = 0.;
	for (i = 0; i < len_vect-2; i++){
		integral +=  (x[i+1] - x[i]) * (vect[i] + vect[i+1]) / 2.;
		if(integral>frac){
			return i;
		}
	}
	return 0;
}

/**
 * \brief Gets a percentil location integrating from the right
 * \param x a pointer to double, the x axis values
 * \param vect a pointer to double, the y(x) axis values
 * \param len_y an int, the length of the vect array
 * \param frac a double, the given percentil
 * \return a double, the value of the percentil location
 *
 * Gets a percentil location integrating from the right
 */
int GetLocation_IntTrap_Right(double* x, double* vect, int len_vect, double frac){
	double integral;
	int i;

	integral = 0.;
	for (i = len_vect-2; i >= 0; i--){
		integral +=  (x[i+1] - x[i]) * (vect[i] + vect[i+1]) / 2.;
		if(integral>frac){
			return i;
		}
	}
	return 0;
}

