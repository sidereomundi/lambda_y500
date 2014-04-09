/**
 * \file MassConversion.h
 * \brief MassConversion module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module for conversion between mass definitions
 */

#ifndef MASSCONVERSION_H_
#define MASSCONVERSION_H_

#include "Function.h"

double MassConversion(Function_2D* f, double, double);
void MassConversion_optimal1D(double*, double*, int, double*, double*, int);
void MassConversion_optimal2D(Function_2D*, double*, int, double, double*);

double InvMassConversion(Function_2D* f, double, double, double, double);
double InvMassConversion_optimal1D(double*, double*, int, double, double, double);
void InvMassConversion_optimal2D(Function_2D*, double*, int, double, double, double*);

#endif /* MASSCONVERSION_H_ */
