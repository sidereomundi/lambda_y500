/**
 * \file ScalingRelation.h
 * \brief ScalingRelation module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module for scaling relations
 */

#ifndef SCALINGRELATION_H_
#define SCALINGRELATION_H_

#include "Function.h"
#include "Parameters.h"                                                                                
#include <SebCosmolib.h>


double Get_M200_M_DELTA(double, double, double);
/**
 * \def SZ_A
 * \brief prior center for SZ normalization
 * \def SZ_B
 * \brief prior center for SZ slope
 * \def SZ_C
 * \brief prior center for SZ redshift evolution
 * \def SZ_D
 * \brief prior center for SZ scatter
 * \def SZ_sigmaA
 * \brief uncertainty on the prior center for SZ normalization
 * \def SZ_sigmaB
 * \brief uncertainty on the prior center for SZ slope
 * \def SZ_sigmaC
 * \brief uncertainty on the prior center for SZ redshift evolution
 * \def SZ_sigmaD
 * \brief uncertainty on the prior center for SZ scatter
 * \def S_A
 * \brief prior center for dispersion normalization
 * \def S_B
 * \brief prior center for dispersion slope
 * \def S_C
 * \brief prior center for dispersion redshift evolution
 * \def S_D
 * \brief prior center for dispersion scatter
 * \def S_sigmaA
 * \brief uncertainty on the prior center for dispersion normalization
 * \def S_sigmaB
 * \brief uncertainty on the prior center for dispersion slope
 * \def S_sigmaC
 * \brief uncertainty on the prior center for dispersion redshift evolution
 * \def S_sigmaD
 * \brief uncertainty on the prior center for dispersion scatter
 * \def X_A
 * \brief prior center for X-ray normalization
 * \def X_B
 * \brief prior center for X-ray slope
 * \def X_C
 * \brief prior center for X-ray redshift evolution
 * \def X_D
 * \brief prior center for X-ray scatter
 * \def X_sigmaA
 * \brief uncertainty on the prior center for X-ray normalization
 * \def X_sigmaB
 * \brief uncertainty on the prior center for X-ray slope
 * \def X_sigmaC
 * \brief uncertainty on the prior center for X-ray redshift evolution
 * \def X_sigmaD
 * \brief uncertainty on the prior center for X-ray scatter
 * \def Zeta_Threshold
 * \brief threshold in zeta below which the scaling relation is modified
 * \def Xi_Threshold
 * \brief threshold in xi below which the scaling relation is modified
 * \def ZetaXi_Slope
 * \brief slope between zeta and xi below the threshold
 */


/****** SZE Scaling relation ******/
#define SZ_A 5.58
#define SZ_B 1.32
#define SZ_C 0.87
#define SZ_D 0.24
/* absolute uncertainties */
#define SZ_sigmaA 1.67
#define SZ_sigmaB 0.26
#define SZ_sigmaC 0.44
#define SZ_sigmaD 0.16

#define Zeta_Threshold 2.6457513110645907
#define Xi_Threshold 3.1622776601683795
#define ZetaXi_Slope 0.8366600265340756

/****** Velocity Dispersion Scaling relation ******/
/*#define S_A 1047.61
#define S_B 0.342
#define S_C 0.319
#define S_D 0.266*/
/* absolute uncertainties */
#define S_sigmaA 53.0
#define S_sigmaB 0.01
#define S_sigmaC 0.02
#define S_sigmaD 0.08

/****** X-ray Scaling relation ******/
#define X_A 5.77
#define X_B 0.57
#define X_C -0.4
#define X_D 0.12
/* absolute uncertainties */
#define X_sigmaA 0.56 /*  9% absolute mass calibration systematics + stat */
#define X_sigmaB 0.03
#define X_sigmaC 0.2 /* 50% uncertainty on the redshift evolution parameter following B11 */
#define X_sigmaD 0.08 /* 66% uncertainty on the scatter following B11 */



/****** SZE Scaling relation ******/
double R500_in_arcmin_given_M500(double M500, double z,                               
				 double h0, double Om, double Ol, double w0, double wa);


Function_1D* Get_P_Y500cyl_obs_lambda(Parameters_struct* params,
				      Function_1D* vect_p_m_lambda_debiased,
				      double redshift, 
				      double* y0,
				      double* dy0,
				      double* filter,
				      int NY500,
				      y0_to_Y500cyl_structure* y0_to_Y500cyl);


Function_1D* Get_P_Y500cyl_lambda(Parameters_struct* params,
				  Function_1D* p_m_lambda, 
				  double redshift);

double m500c2Y500_spher(double m500c, 
			double z,                                                                                    
                        double h0, 
			double Om, 
			double Ol, 
			double w0, 
			double wa);
double m500c2Y500_cyl(double m500c, 
			double z,                                                                                    
                        double h0, 
			double Om, 
			double Ol, 
			double w0, 
			double wa);
double Xi2Zeta(double);
double Zeta2Xi(double);
double Zeta2Mass(double zeta, 
		 double z, 
		 double Asz, 
		 double Bsz, 
		 double Csz, 
		 double h0, 
		 double Om, 
		 double Ol, 
		 double w0, 
		 double wa);
double Mass2Zeta(double mass, 
		 double z, 
		 double Asz, 
		 double Bsz, 
		 double Csz, 
		 double h0, 
		 double Om, 
		 double Ol, 
		 double w0, 
		 double wa);

double** Get_array_lambda_M(double* lambda_vect, int len_lambda_vect, double Dlambda);


double** Get_Plan_xi_M(double* xi_vect, 
		       int len_xi_vect, 
		       double Dsz);

Function_1D* Get_P_xi_M(double** plan_xi_m, 
			double* xi_vect, 
			int len_xi_vect, 
			double xi, 
			double redshift);

Function_1D* Get_P_M_lambda(Parameters_struct* params,
			    Function_1D* P_lambda0_lambdaobs,
			    Function_2D* MassFunc,
			    double redshift);

double lambda2mass(double lambda, 
		   double z,double Alambda, 
		   double Blambda, double Clambda,
		   double h0, double Om, 
		   double Ol, double w0, double wa);


Function_1D* Get_P_lambda0_lambdaobs(double** plan_lambda_m, 
				     double* lambda_vect, 
				     int len_lambda_vect, 
				     double lambda, double lambdaerr);

Function_1D* Get_P_M_xi(Function_1D* P_xi_M, 
			Function_2D* MassFunc, 
			double redshift, 
			double Asz, 
			double Bsz, 
			double Csz,
			double h0, 
			double Om, 
			double Ol, 
			double w0, 
			double wa);

/****** Velocity Dispersion Scaling relation ******/
double Disp2Mass(double, 
		 double, 
		 double, 
		 double,
		 double,
		 double,
		 double,
		 double,
		 double,
		 double);
double Mass2Disp(double, 
		 double, 
		 double, 
		 double,
		 double,
		 double,
		 double,
		 double,
		 double,
		 double);
double Disp_IntScat(double z, double Ds);
Function_1D* Get_P_Sigma_xi(Function_1D* p_m_xi,
			    Function_2D* m200c_2_m500c,
			    double redshift, 
			    double veldisp,
			    double veldisperr, 
                int Ngal,
			    double* veldisp_vect, 
			    double As,
			    double Bs, 
			    double Cs,
			    double Ds,
                double Ds_0,
                double Ds_N,
			    double h0, 
			    double Om, 
			    double Ol, 
			    double w0, 
			    double wa);
Function_1D* Get_P_Sigma_Yx(Function_1D* p_m_yx,
			    double redshift, 
			    double veldisp, 
			    double veldisperr, 
			    double As,
			    double Bs, 
			    double Cs,
			    double Ds,
			    double h0, 
			    double Om, 
			    double Ol, 
			    double w0, 
			    double wa);

/****** X-ray Scaling relation ******/
double Yx2Mass(double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double);
double Mass2Yx(double,
	       double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double, 
	       double,
	       double,
	       double);
Function_1D* Get_P_Yx_xi(Function_1D*, 
			 double, 
			 double, 
			 double,
			 double*,
			 double, 
			 double, 
			 double,
			 double,
			 double, 
			 double, 
			 double, 
			 double, 
			 double,
			 short);
Function_1D* Get_P_M200c_Yx(Function_2D*,
			    double, 
			    double, 
			    double, 
			    double*,
			    int,
			    double, 
			    double,
			    double,
			    double,
			    double,
			    double,
			    double,
			    double,
			    double);




#endif /* SCALINGRELATION_H_ */

