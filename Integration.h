/**
 * \file Integration.h
 * \brief Integration module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to integrate a function
 */

#ifndef INTEGRATION_H_
#define INTEGRATION_H_

double IntTrap_conststep(double, double*, int);
double IntTrap_varstep(double*, double*, int, int);

int GetLocation_IntTrap_Left(double*, double*, int, double);
int GetLocation_IntTrap_Right(double*, double*, int, double);

#endif /* INTEGRATION_H_ */
