/**
 * \file Statistics.h
 * \brief Statistics module header
 * \date on: 11 August 2011
 * \author Gurvan
 *
 * Module to calculate the likelihood
 */

#include "Function.h"

#ifndef STATISTICS_H_
#define STATISTICS_H_

/**
 * \def EPS1
 * \brief epsilon for convergence
 * \def EPS2
 * \brief epsilon for convergence
 * \def WRONG_PULL
 * \brief wrong value for the pull
 */
#define EPS1 0.001
#define EPS2 1.e-8
#define WRONG_PULL -1000

/**
 * \struct erf_struct
 * \brief Error function structure
 */
typedef struct {
	long double* x;
	long double* y;
} erf_struct;


erf_struct* InitErf(char*);
void FreeErf(erf_struct**);
double inverf(erf_struct*, long double);
double Get_Pull(Function_1D*, double);
double Get_GaussPull(Function_1D*, double, erf_struct*);
double Get_LogPull(Function_1D*, double);
double Get_KS_p(double);
double KS_Pvalue(double*, int);


#endif /* STATISTICS_H_ */
