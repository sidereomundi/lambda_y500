#ifndef _COSMOROUTINES_H_
#define _COSMOROUTINES_H_

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define Squ(x) ((x)*(x))
#define Cube(x) ((x)*(x)*(x))
#define Quad(x) ((x)*(x)*(x)*(x))

#ifndef PI
#define PI 3.14159265359
#endif

// Physics Constants
#ifndef DIST_H
#define DIST_H 2997.92458   /* c/H*h in Mpc*/
#endif

struct dAparams {
    double Om,Ol,w0,wa;
};


double SebCosmolib_E_z(double, double, double, double);
double SebCosmolib_integrand_Ez(double, void *);
void SebCosmolib_CorrectCurvature(double, double *);
double SebCosmolib_AngDiamDist(double, double, double, double);


#endif
