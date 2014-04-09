/**
 * \file Cosmology.h
 * \brief Cosmology module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to calculate cosmological quantities and distances, also contains functions for the calculation of the Tinker mass function
 */

#ifndef COSMOLOGY_H_
#define COSMOLOGY_H_

/** 
 * \def c 
 * \brief the speed of light in km s^-1
 * \def H0
 * \brief the Hubble constant in km s^-1 Mpc^-1
 * \def H0_Mpc
 * \brief the Hubble constant in Mpc s^-1 Mpc^-1
 * \def h_o
 * \brief H0/100
 * \def OmegaM
 * \brief Matter density
 * \def OmegaL
 * \brief Dark Energy density
 * \def G
 * \brief the gravitational constant in m^3 kg^-1 s^-2
 * \def G_Mpc_Msun
 * \brief the gravitational constant in Mpc^3 Msun^-1 s^-2
 * \def W0
 * \brief 0 order parameter in the dark energy equation of state
 * \def WA
 * \brief 1st order parameter in the dark energy equation of state
 * \def sigmah_0
 * \brief uncertainty of h_0
 * \def sigmaOmegaM
 * \brief uncertainty of OmegaM
 * \def Mpc_2_m
 * \brief conversion factor between Mpc and m
 * \def Msun_2_kg
 * \brief conversion factor between Msun and kg
 */

//#define c 299792.458
#define H0 71.7
#define H0_Mpc 2.2750250937374343e-18
#define h_0 0.717
#define OmegaM 0.255
#define OmegaL 0.745
#define G_Mpc_Msun 4.5173701455761134e-48

#ifndef DIST_H
#define DIST_H 2997.92458   /* c/H*h in Mpc*/
#endif
#ifndef PI
#define PI 3.14159265359
#endif
#ifndef G
#define G 6.67384e-11       /* Gravitational constant in m^3 kg^-1 s^-2*/
#endif
#ifndef Mpc_in_m
#define Mpc_in_m 3.08567758e+22  /* conversion from Mpc to m */
#endif
#ifndef MSOLAR
#define MSOLAR 1.9891e+30   /* Solar mass in kg */
#endif

#define W0 -1.0
#define WA 0.0

#define sigmah_0 0.016
#define sigmaOmegaM 0.016

#define Mpc_2_m 3.08567758e22 /* m */
#define Msun_2_kg 1.9891e30 /* kg */

/*
double h(double);
double h_z(double, double, double, double, double, double);
double E_z_masscalib(double, double, double, double, double);
double E(double);
double integrand_Ez_masscalib(double,void *);
double dA(double,double, double,double,double,double);
double dL(double);
*/
double W(double);
double dWdx(double);
double R(double);
double dRdM(double);
double f_sigma(double, double, double);
double sigma2(double*, double*, int, double, double);
double dsigma2dM(double*, double*, int, double, double);
double dndlnM(double*, double*, int, double, double, double);

#endif /* COSMOLOGY_H_ */
