/**
 * \file Interpolation.h
 * \brief Interpolation module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to interpolate functions
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

double Interp1D_Scalar(double*, double*, int, double);
void Interp1D_Vector(double*, double*, int, double*, double*, int);

double Interp2D_Scalar(double*, double*, double**, int, int, double, double);

#endif /* INTERPOLATION_H_ */
