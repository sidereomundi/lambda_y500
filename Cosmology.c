/**
 * \file Cosmology.c
 * \brief Cosmology module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to calculate cosmological quantities and distances, also contains functions for the calculation of the Tinker mass function
 */
#include "Cosmology.h"
#include "Integration.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_integration.h>


/**
 * \brief Calculates h(z)
 * \param z0 a double, the redshift
 * \return h(z) a double
 *
 * Calculates h(z) for a fixed LCDM cosmology defined by h_0, OmegaM and OmegaL in the code
 */
double h(double z0){
	return h_0  * sqrt(OmegaM * pow(1+z0, 3) + OmegaL);
}

/**
 * \brief Calculates h(z)
 * \param z a double, the redshift
 * \param h0 a double, the Hubble constante divided by 100
 * \param Om a double, OmegaM
 * \param Ol a double, OmegaL
 * \param w0 a double, 0 order parameter in the dark energy equation of state
 * \param wa a double, 1st order parameter in the dark energy equation of state
 * \return h(z) a double
 *
 * Calculates h(z) for a given cosmology
 */
double h_z(double z, double h0, double Om, double Ol, double w0, double wa){
  return h0 * E_z_masscalib(z, Om, Ol, w0, wa);
}


/**
 * \brief Calculates E(z)
 * \param z0 a double, the redshift
 * \return E(z) a double
 *
 * Calculates E(z) for a fixed LCDM cosmology defined by h_0, OmegaM and OmegaL in the code
 */
double E(double z0){
	return h(z0)/h_0;
}

/**
 * \brief Calculates E(z)
 * \param z a double, the redshift
 * \param Om a double, OmegaM
 * \param Ol a double, OmegaL
 * \param w0 a double, 0 order parameter in the dark energy equation of state
 * \param wa a double, 1st order parameter in the dark energy equation of state
 * \return E(z) a double
 *
 * Calculates E(z) for a given cosmology
 */
double E_z_masscalib(double z, double Om, double Ol, double w0, double wa){
    double integ_w = (1.+w0) * log(1.+z);// - z/(1.+z)*wa;
  return sqrt(Om*pow(1.+z, 3.) + (1.-Om-Ol)*(1.+z)*(1.+z) + Ol*exp(3.*integ_w));
}

// Curvature correction for non-flat cosmology
void correctCurvature_masscalib(double omegaK, double *x)
{
	if (omegaK>1.e-4)
		*x *= sinh(sqrt(omegaK))/sqrt(omegaK);
	else if (omegaK<-1.e-4)
		*x *= sin(sqrt(-omegaK))/sqrt(-omegaK);
}


// returns the angular diameter distance
double dA(double z, double h, double Om, double Ol, double w0, double wa)
{
    struct dAparams dApar;
    double x,error;
    
    dApar.Om = Om;
    dApar.Ol = Ol;
    dApar.w0 = w0;
    dApar.wa = wa;
    
    // GSL integration init
    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);
    gsl_function F;
    F.function = &integrand_Ez_masscalib;
    F.params = &dApar;
    
    gsl_integration_qag (&F, 0., z, 0, 1.e-6, 1000, 6, w, &x, &error);
    gsl_integration_workspace_free (w);
    
    // los comoving distance
    x *= DIST_H/h;
    
    // transverse comoving distance
	correctCurvature_masscalib(1. -Om -Ol, &x);
    
    // angular diameter distance
	return x/(1.+z);
}




// integrand for determining the comoving distance
double integrand_Ez_masscalib(double z, void * params) {
    struct dAparams *p = (struct dAparams*) params;
	return 1./E_z_masscalib(z,p->Om,p->Ol,p->w0,p->wa);
}
/**
 * \brief Calculates dL(z)
 * \param z0 a double, the redshift
 * \return dL(z) a double
 *
 * Calculates the luminosity distance for a fixed LCDM cosmology defined by h_0, OmegaM and OmegaL in the code
 */
//double dL(double z0){
//	return dA(z0) * (1.+z0) * (1.+z0);
//}




/* for mass function */
double rho_mean(){ /* Msun Mpc-3 */
	return 3.*H0_Mpc*H0_Mpc*OmegaM / (8.*M_PI*G_Mpc_Msun);
}

double W(double x){
	if (x < 5.e-2) return 1. - x*x/10. + x*x*x*x/280. - x*x*x*x*x*x /134400.;
	return (3./pow(x, 3)) * (sin(x) - x*cos(x));
}

double dWdx(double x){
	if (x < 3.e-2) return -x/5. + x*x*x/70. + x*x*x*x*x/22400.;
	return (3./pow(x, 4)) * (3.*x*cos(x) + (x*x-3.)*sin(x));
}

double R(double M){ /* Mpc, M in Msun */
	return pow(M/(4.*M_PI*rho_mean()/3.), 1./3.);
}

double dRdM(double M){
	return R(M)/(3.*M);
}

double f_sigma(double sigma, double z, double delta){
	double A, A0;
	double a0, b0;
	double a, b, c0;
	double alpha;
	double logdelta;

	if(delta>=1600)
		A0 = 0.26;
	else
		A0 = 0.1*log10(delta) - 0.05;

	alpha = pow(10., -pow(0.75/log10(delta/75.), 1.2));

	logdelta = log10(delta);
	a0 = 1.43 + pow(logdelta-2.3, 1.5);
	b0 = 1.00 + pow(logdelta-1.6, -1.5);
	c0  = (logdelta>2.35 ? 1.20 + pow(logdelta-2.35, 1.6): 1.19);

	A = A0 * pow(1.+z, -0.14);
	a = a0 * pow(1.+z, -0.06);
	b = b0 * pow(1.+z, -alpha);

	return A * (pow(sigma/b, -a) + 1) * exp(-c0/(sigma*sigma));
}

double sigma2(double* vect_lnk, double* vect_Pk, int len_lnk, double z, double M){
	double* integrand;
	double k;
	int i;
	double WkR;
	double RM;
	double integral;

	RM = R(M);

	/* allocate memory for the integration */
	integrand = (double*) malloc(len_lnk * sizeof(double));

	for(i = 0; i < len_lnk; i++){
		k = exp(vect_lnk[i]);
		WkR = W(k*RM);
		integrand[i] = vect_Pk[i] * WkR*WkR * k*k*k;
	}
	integral = IntTrap_varstep(vect_lnk, integrand, len_lnk, len_lnk);

	/* free */
	free(integrand);

	return integral / (2.*M_PI*M_PI);
}

double dsigma2dM(double* vect_lnk, double* vect_Pk, int len_lnk, double z, double M){
	double* integrand;
	double k;
	int i;
	double RM;
	double dRM;
	double kR;
	double integral;

	RM = R(M);
	dRM = dRdM(M);

	/* allocate memory for the integration */
	integrand = (double*) malloc(len_lnk * sizeof(double));

	for(i = 0; i < len_lnk; i++){
		k = exp(vect_lnk[i]);
		kR = k*RM;
		integrand[i] = vect_Pk[i] * W(kR) * dWdx(kR) * dRM * k*k*k*k;
	}
	integral = IntTrap_varstep(vect_lnk, integrand, len_lnk, len_lnk);

	/* free */
	free(integrand);

	return 2. * integral / (2.*M_PI*M_PI);

}

double dndlnM(double* vect_lnk, double* vect_Pk, int len_lnk, double M, double z, double delta){
	double sigma, sigmasq;
	double dsigma2dlnM;
	double result;

	sigmasq = sigma2(vect_lnk, vect_Pk, len_lnk, z, M);
	sigma = sqrt(sigmasq);
	dsigma2dlnM = M * dsigma2dM(vect_lnk, vect_Pk, len_lnk, z, M);

	result = -f_sigma(sigma, z, delta) * rho_mean() * dsigma2dlnM / (M * 2. * sigmasq);

	return result;

}
