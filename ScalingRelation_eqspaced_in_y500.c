/**
 * \file ScalingRelation.c
 * \brief ScalingRelation module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module for scaling relations
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "SZFitC.h"
#include "utils.h"
#include "Data.h"
#include "Parameters.h"
#include "ScalingRelation.h"
#include "Cosmology.h"
#include "Convolution.h"
#include "Interpolation.h"
#include "Integration.h"
#include "MassConversion.h"

#include <SebCosmolib.h>

// Mass observable relations

double lambda2mass(double lambda, double z,
		double Alambda, double Blambda, double Clambda,
		   double h0, double Om, double Ol, double w0, double wa){
  //	return pow(lambda/(Alambda), 1./Blambda);
  return pow(lambda/(Alambda * pow(SebCosmolib_E_z(z, Om, Ol, w0)/ SebCosmolib_E_z(0.5, Om, Ol, w0),Clambda)), 1./Blambda);
}


double Xi2Zeta(double xi){
  double zeta;
/*  if (xi > Xi_Threshold){
    zeta = sqrt(xi*xi-3.);
  }else{
    zeta = ZetaXi_Slope * xi;
  }*/
    // Adapted crazy SPT relation according to Tijmen
    if (xi>sqrt(7.))
        zeta = sqrt(xi*xi-3.);
    else
		zeta = sqrt(4./7.)*xi;
  return zeta;
}


double Zeta2Xi(double zeta){
  double xi;
/*  if (zeta > Zeta_Threshold){
    xi = sqrt(zeta*zeta+3);
  }else{
    xi = zeta / ZetaXi_Slope;
  }*/
    // Crazy SPT relation according to Tijmen
    if(zeta>2.)
        xi = sqrt(zeta*zeta+3.);
    else
        xi = zeta;
  return xi;
}


double Zeta2Mass(double zeta, double z, 
		double Asz, double Bsz, double Csz, double h0, 
		double Om, double Ol, double w0, double wa){
	return (3./h0) * pow(zeta/(Asz*pow(SebCosmolib_E_z(z, Om, Ol, w0)/SebCosmolib_E_z(0.6, Om, Ol, w0), Csz)), (1./Bsz));
}


double Mass2Zeta(double mass, double z, 
		double Asz, double Bsz, double Csz, 
		double h0, double Om, double Ol, double w0, double wa){
	return Asz*pow(mass*h0/3., Bsz) * pow(SebCosmolib_E_z(z, Om, Ol, w0)/SebCosmolib_E_z(0.6, Om, Ol, w0), Csz);
}

// Arnaud et al. 2010 
double m500c2Y500_spher(double m500c, double z, 
			double h0, double Om, double Ol, double w0, double wa){
  double h70=h0/70.;
  return pow(10,-4.739)*pow(m500c/3.e14*h70,1.790)*pow(h70,-5./2)*pow(h70*SebCosmolib_E_z(z, Om, Ol, w0),2./3);
}

double m500c2Y500_cyl(double m500c, double z, 
			double h0, double Om, double Ol, double w0, double wa){
  return 1.22*m500c2Y500_spher(m500c,z,h0,Om,Ol,w0,wa);
}


double Yx2Mass(double yx, double z, 
		double Ax, double Bx, double Cx, 
		double h0, double Om, double Ol, double w0, double wa){
	return Ax * sqrt(h0) * pow(yx/3., Bx) * pow(SebCosmolib_E_z(z, Om, Ol, w0), Cx);
}


double Mass2Yx(double mass, double z, 
		double Ax, double Bx, double Cx,
		double h0, double Om, double Ol, double w0, double wa){
	return 3. * pow(mass / (Ax * sqrt(h0) * pow(SebCosmolib_E_z(z, Om, Ol, w0), Cx)) , 1./Bx);
}


///////////////////////////////


/* Provide the array containing the relation between lambda_0 and lambda_obs
	For a given underlying lambda_0, apply the scatter to provide P(lambda_obs|lambda_0)*/
double** Get_array_lambda_M(double* lambda_vect, int len_lambda_vect, double Dlambda)
{
	int i,j;
	double lambda_0;
	
	// Allocate memory for the output, shape is len_lambda_vect*len_lambda_vect
	double** result = (double**) malloc(len_lambda_vect*sizeof(double*));
	for(i=0;i<len_lambda_vect;i++)
		result[i] = (double*) malloc(len_lambda_vect*sizeof(double));
	
	// For each lambda, apply the lognormal uncertainty
	// definition: lognormal(x, mu, sigma)
	for(j = 0; j < len_lambda_vect; j++)
	{
		lambda_0 = lambda_vect[j];
		for(i = 0; i < len_lambda_vect; i++)
			result[j][i] = lognormal(lambda_vect[i], lambda_0, Dlambda);
	}
	
	// Each result[j][i] contains P(lambda_obs[i]|lambda_0[j])
	
	return result;
}




/* Return P(M|lambda) (not yet debiased)
	Basically, you read plan_lambda_m in the "wrong" direction:
	For a given lambda_observed, read in the lambda_0 direction*/
Function_1D* Get_P_lambda0_lambdaobs(double** plan_lambda_m, 
			double* lambda_vect, 
			int len_lambda_vect, 
			double lambda){
  int i, k1, k2;

  // Allocate memory for the output
  Function_1D* result = (Function_1D*) malloc(sizeof(Function_1D));
  result->N = len_lambda_vect;
  result->x = (double*) malloc(len_lambda_vect * sizeof (double));
  result->y = (double*) malloc(len_lambda_vect * sizeof (double));

  // find where lambda is in lambda_vect
  if(lambda < lambda_vect[0] || lambda > lambda_vect[len_lambda_vect-1]){
    printf("Get_P_xi_M::error, xi= %lf is out of range: %lf - %lf\n", lambda, lambda_vect[0], lambda_vect[len_lambda_vect-1]);
  }
  i = 0;
  while(1){
    if(lambda < lambda_vect[i] || i > len_lambda_vect-2) break;
    i++;
  }
  k1 = i-1;
  k2 = i;
  
  // interpolate P(lambda|M)
  for(i = 0; i < len_lambda_vect; i++)
  {
    result->x[i] = lambda_vect[i];
    result->y[i] = plan_lambda_m[i][k1] + (lambda-lambda_vect[k1]) * (plan_lambda_m[i][k2]-plan_lambda_m[i][k1])/(lambda_vect[k2]-lambda_vect[k1]);
  }

  return result;
}

// Multiply P(M|lambda)_0 with the mass function and return P(M|lambda)_debiased
Function_1D* Get_P_M_lambda (Parameters_struct* params,
			    Function_1D* P_lambda0_lambdaobs, 
			    Function_2D* MassFunc, 
			    double redshift){
	int n = P_lambda0_lambdaobs->N;
	double* m500c;
	double* dNdM;
	double* dNdM_redshift;
	int i;
	FILE* fp;
	
	// Get parameter values
	double Alambda = params->theta[I_lambda_A];
	double Blambda = params->theta[I_lambda_B];
	double Clambda = params->theta[I_lambda_C];
	double h0 = params->theta[I_H0];
	double Om = params->theta[I_OM];
	double Ol = params->theta[I_OL];
	double w0 = params->theta[I_W0];
	double wa = params->theta[I_WA];

	// Allocate memory
	dNdM = malloc(n*sizeof(double));
	dNdM_redshift = malloc(MassFunc->Nx*sizeof(double));

	// Allocate the returned array, shape is n*n
	Function_1D* result = (Function_1D*) malloc(sizeof(Function_1D));
	result->N = n;
	result->x = (double*) malloc(n * sizeof (double));
	result->y = (double*) malloc(n * sizeof (double));

	// Interpolate the MF to the cluster's redshift
	for(i = 0; i<MassFunc->Nx; i++)
    	dNdM_redshift[i] = Interp1D_Scalar(MassFunc->y, MassFunc->z[i], MassFunc->Ny, redshift);

	//	for(i = 0; i<MassFunc->Nx; i++)printf("%e\n",MassFunc->z[i][0]);
	//for(i = 0; i<MassFunc->Nx; i++)printf("%e\n",dNdM_redshift[i]);
	//return 0;

    
	// Convert lambda to mass
	for(i = 0; i < n; i++)
		result->x[i] = lambda2mass(P_lambda0_lambdaobs->x[i], redshift, Alambda, Blambda, Clambda, h0, Om, Ol, w0, wa);
	//	printf("%f %f %f %f %f %f %f %f %f \n",redshift, Alambda, Blambda, Clambda, h0, Om, Ol, w0, wa);
	//	for(i = 0; i < n; i++)printf("%f %f %e\n",P_lambda0_lambdaobs->x[i],result->x[i],result->y[i]);
	//	return 0;


	// Interpolate the mass function (correct redshift) to each of the masses just calculated
	Interp1D_Vector(MassFunc->x, dNdM_redshift, MassFunc->Nx, result->x, dNdM, n);

	// Now do it: P(M|lambda)*P(M)
	for(i = 0; i < n; i++)
		result->y[i] = P_lambda0_lambdaobs->y[i] * dNdM[i];

    free(dNdM);
    free(dNdM_redshift);

	return result;
}



// Calculate P(zeta|lambda)
Function_1D* Get_P_zeta_lambda
(Parameters_struct* params,
 Function_1D* p_m_lambda, 
 double redshift, 
 double zeta){
  // Get parameter values
	double Asz = params->theta[I_SZ_A];
	double Bsz = params->theta[I_SZ_B];
	double Csz = params->theta[I_SZ_C];
	double Dsz = params->theta[I_SZ_D];
	double h0 = params->theta[I_H0];
	double Om = params->theta[I_OM];
	double Ol = params->theta[I_OL];
	double w0 = params->theta[I_W0];
	double wa = params->theta[I_WA];

	int i, j;
	int n = p_m_lambda->N;
	double zeta_vect[n];
	double tmp[n];
	double dzeta;
	double maxvalue;
	Function_1D freg;
	double dzetadm;
	double* lognormal_zeta = (double*) malloc(n * sizeof(double));
	double* dpdzeta = (double*) malloc(n * sizeof(double));
	double* dpdzeta_eqspaced = (double*) malloc(n * sizeof(double));
	double* normal_zeta = (double*) malloc(n * sizeof(double));
	double* zeta_eqspaced = (double*) malloc(n * sizeof(double));
	double* zeta_distrib = (double*) malloc(n * sizeof(double));
	double* integrand = (double*) malloc(n * sizeof(double));
	double** array_2d = (double**) malloc(n * sizeof(double*));
	Function_1D* fout = (Function_1D*) malloc(sizeof(Function_1D));
	for(i = 0; i < n; i++) array_2d[i] = (double*) malloc(n * sizeof(double));
	fout->x = (double*) malloc(n * sizeof(double));
	fout->y = (double*) malloc(n * sizeof(double));

  	// zeta(M)
    for (i = 0 ; i < n; i++)
        zeta_vect[i] = Mass2Zeta(p_m_lambda->x[i], redshift, Asz, Bsz, Csz, h0, Om, Ol, w0, wa);
    
	// Calculate dpdzeta from dpdm
	for(i = 0 ; i < n-1; i++){
		if (i < n-1) dzetadm = (zeta_vect[i+1]-zeta_vect[i])/(p_m_lambda->x[i+1]-p_m_lambda->x[i]);
		dpdzeta[i] = p_m_lambda->y[i] / dzetadm;
	}
	dpdzeta[n-1] = 0.;

	// equally-spaced zeta grid
	double zeta_max = minimum(zeta_vect[n-1],abs(zeta*10.));
	if(zeta_max<10)zeta_max = 10.;
	if(zeta_max<zeta+1)zeta_max = zeta+1;
	double zeta_min = -5+1.e-3;
	dzeta = (zeta_max - zeta_min) / ((double) n);

	for(i = 0 ; i < n; i++)
	  zeta_eqspaced[i] = zeta_min+((double)i)*dzeta;

/* 	if(zeta>zeta_eqspaced[n-1]){ */
/* 	  printf("WARNING: zeta %lf > zeta_eqspaced %lf\n",zeta,zeta_eqspaced[n-1],p_m_lambda->x[n-1]); */
/* 	  printf("WARNING: zeta min %lf, zeta_max %lf\n",zeta_min,zeta_max); */
/* 	  //	  for(i=0;i<n;i++)printf("%lf %lf %lf\n",zeta_eqspaced[i],zeta_vect[i],p_m_lambda->x[i]); */
/* 	} */
/*     if(zeta<zeta_eqspaced[0]) */
/*       printf("WARNING: zeta %lf < zeta_eqspaced %lf\n",zeta,zeta_eqspaced[0],p_m_lambda->x[0]); */
    
	// Interpolate dpdzeta to match the values in zeta_eqspaced and put it into dpdzeta_eqspaced
    Interp1D_Vector(zeta_vect, dpdzeta, n, zeta_eqspaced, dpdzeta_eqspaced, n);
    //AS *** set to zero for neg zeta or to 1? Shouldn't matter?
    for(i=0;i<n;i++)if(zeta_eqspaced[i]<0)dpdzeta_eqspaced[i]=0.;
    
    /* Account for the SZ scatter
       convolve a lognormal distribution with a normal distribution*/
    double zeta_0;
    for (i = 0 ; i < n; i++)
      {
	zeta_0 = zeta_eqspaced[i];
	
        for (j = 0 ; j < n; j++)
	  {
	    lognormal_zeta[j] = 0;
	    if(zeta_0>0)
	      if(zeta_eqspaced[j]>0)
		lognormal_zeta[j] = lognormal(zeta_eqspaced[j], zeta_0, Dsz);
	    // AS Which one???
	    //	    normal_zeta[j] = normal(zeta_eqspaced[j], 0., 1.);
	    	    normal_zeta[j] = normal(zeta_eqspaced[j], zeta_eqspaced[n/2], 1.);
	  }
	
        // fix the zero problem
        maxvalue = 1.e-10;
        for(j = 0; j < n; j++){
            if(lognormal_zeta[j] < maxvalue) lognormal_zeta[j] = maxvalue;
			if(normal_zeta[j] < maxvalue) normal_zeta[j] = maxvalue;
        }

        // get the zeta distribution for each initial mass
	//	if(i==n/2)for(j=0;j<n;j++)printf("%d %f %f %f \n",i,zeta_eqspaced[j],normal_zeta[j], lognormal_zeta[j]);		    

	Convolve(normal_zeta, lognormal_zeta, n, n, tmp);
  
	   // Regularization
       freg.N = n;
       freg.x = zeta_eqspaced;
       freg.y = tmp;
	// AS remove regularization!!!! ******
	//        Regularize_1D(&freg, zeta_vect[i], REGUL_LOG);
        
        // copy in 2D array
		// array_2d[mu][x]
       for (j = 0 ; j < n; j++)
	 array_2d[i][j] = freg.y[j];
      }
    
    // Calculate the zeta distribution
    for (i = 0 ; i < n; i++)
      {
        for (j = 0 ; j < n; j++)
	  {
	    integrand[j] = array_2d[j][i] * dpdzeta_eqspaced[j];	    
	  }
	//	printf("%f %f \n",integrand[i],dpdzeta_eqspaced[i]);	

        zeta_distrib[i] = IntTrap_varstep(zeta_eqspaced, integrand, n, n);
      }
    
    fout->N = n;
    for (i = 0 ; i < n; i++) {
        fout->x[i] = zeta_eqspaced[i];
        fout->y[i] = zeta_distrib[i];
	//	if(zeta<1.115 && zeta>1.1147)
	//	printf("%f %f %f %f \n",p_m_lambda->x[i],p_m_lambda->y[i],zeta_eqspaced[i],zeta_distrib[i]);
    }
    

	free(integrand);
	free(lognormal_zeta);
	free(dpdzeta);
	free(dpdzeta_eqspaced);
	free(normal_zeta);
	free(zeta_eqspaced);
	for(i = 0; i < n; i++) free(array_2d[i]);
	free(array_2d);
	free(zeta_distrib);

	return fout;
}


double R500_in_arcmin_given_M500(double M500, double z,                               
				 double h0, double Om, double Ol, double w0, double wa) 
{                                                                                                        
#define RHOCRIT 2.775362e11     // critical density in h^2 M_sol/Mpc^3                                   
  double r500_in_Mpc,angdist,rho_c_z;
  angdist = SebCosmolib_AngDiamDist(z, Om, Ol, w0) / h0; //angular diameter distance in units Mpc/h 
  rho_c_z = RHOCRIT * pow(h0*SebCosmolib_E_z(z,Om,Ol,w0), 2.); //critical density at redshift z:         
  r500_in_Mpc = pow(3.*M500*1.e14/(4.*PI*500.*rho_c_z), 1./3.); // M500 in units 1e14 Msolar                   
  //  printf(" %f %f %f %f %f %f \n",M500,z,angdist,rho_c_z,r500_in_Mpc, r500_in_Mpc/angdist *60.*180./PI);
  return r500_in_Mpc/angdist *60.*180./PI;                                                             
}                                                                                                        


/* Return P(M|lambda) (not yet debiased)
	Basically, you read plan_lambda_m in the "wrong" direction:
	For a given lambda_observed, read in the lambda_0 direction*/
Function_1D* Get_P_Y500cyl_obs_lambda(Parameters_struct* params,
				      Function_1D* vect_p_m_lambda_debiased,
				      double redshift, 
				      double* y0,
				      double* dy0,
				      double* filter,
				      int NY500,
				      y0_to_Y500cyl_structure* y0_to_Y500cyl)
{
  int i, j;
  int Nfilt=params->theta[I_Nfilt];
  
  // Allocate memory for the output
  Function_1D* result = (Function_1D*) malloc(sizeof(Function_1D));
  result->N = NY500;
  result->x = (double*) malloc(NY500 * sizeof (double));
  result->y = (double*) malloc(NY500 * sizeof (double));

  double* R500c_vect = (double*) malloc(NY500 * sizeof(double)); 
  double* y0vect = (double*) malloc(NY500 * sizeof(double));                                                                    
  double* dy0vect = (double*) malloc(NY500 * sizeof(double));                                                                    
  double* y0_to_Y500cyl_factor = (double*) malloc(NY500 * sizeof(double));      
  
  // Allocate memory for the plane NY500 x N500
  double** array_2d = (double**) malloc(NY500*sizeof(double*));
  for(i=0;i<NY500;i++)
    array_2d[i] = (double*) malloc(NY500*sizeof(double));

  double h0 = params->theta[I_H0];
  double Om = params->theta[I_OM];
  double Ol = params->theta[I_OL];
  double w0 = params->theta[I_W0];
  double wa = params->theta[I_WA];
  double Y500cyl_max,Y500cyl_min;
  double central_Y500[NY500];
  double width_of_gaussian[NY500];

  for (i=0;i<NY500;i++)  R500c_vect[i] = R500_in_arcmin_given_M500(vect_p_m_lambda_debiased->x[i],redshift,h0,Om,Ol,w0,wa);


  Interp1D_Vector(filter,y0,Nfilt,R500c_vect,y0vect,NY500);   
  Interp1D_Vector(filter,dy0,Nfilt,R500c_vect,dy0vect,NY500);   
  Interp1D_Vector(y0_to_Y500cyl->r500,y0_to_Y500cyl->y0_to_Y500cyl,146,R500c_vect,y0_to_Y500cyl_factor,NY500);   


  Y500cyl_max = y0vect[0]*y0_to_Y500cyl_factor[0] + 5*dy0vect[0]*y0_to_Y500cyl_factor[0];
  Y500cyl_min = y0vect[0]*y0_to_Y500cyl_factor[0] - 5*dy0vect[0]*y0_to_Y500cyl_factor[0];
  for(j = 0; j < NY500; j++)
    {
      central_Y500[j] = y0vect[j]*y0_to_Y500cyl_factor[j];
    width_of_gaussian[j] = dy0vect[j]*y0_to_Y500cyl_factor[j];
    
    if(central_Y500[j] + 7*width_of_gaussian[j] > Y500cyl_max) 
      Y500cyl_max = central_Y500[j] + 7*central_Y500[j];
    if(central_Y500[j] - 7*width_of_gaussian[j] < Y500cyl_min) 
      Y500cyl_min = central_Y500[j] - 7*central_Y500[j];
    }
  
  for(j = 0; j < NY500; j++)
    {
      result->y[j] = 0;
      for(i = 0; i < NY500; i++)
	{
	  if (j==0) result->x[i] = Y500cyl_min + i*(Y500cyl_max - Y500cyl_min)/(NY500 - 1.);
	  array_2d[j][i] = vect_p_m_lambda_debiased->y[j]*normal(result->x[i],central_Y500[j],width_of_gaussian[j]);
	}
    }
  //

  for(j = 0; j < NY500; j++)
    {
      for(i = 0; i < NY500; i++)
	result->y[j] += array_2d[i][j];
    }

  //  for (i=0;i<NY500;i++){printf(" %f %f \n",result->x[i],result->y[i]);}


  
  free(R500c_vect);
  free(y0vect);
  free(dy0vect);
  free(y0_to_Y500cyl_factor);
  for(i = 0; i < NY500; i++) free(array_2d[i]);
  free(array_2d);
  return result;
}





// Calculate P(zeta|lambda)
Function_1D* Get_P_Y500cyl_lambda(Parameters_struct* params,
				  Function_1D* p_m_lambda, 
				  double redshift){
  // Get parameter values
  double DY500 = params->theta[I_Y500cyl_D];
  double h0 = params->theta[I_H0];
  double Om = params->theta[I_OM];
  double Ol = params->theta[I_OL];
  double w0 = params->theta[I_W0];
  double wa = params->theta[I_WA];

    int i, j;
    int n = p_m_lambda->N;
    double Y500cyl_vect[n];
    double tmp[n];
    double dY500cyl;
    double maxvalue;
    Function_1D freg;
    double dY500cyldm;
    double* lognormal_Y500cyl = (double*) malloc(n * sizeof(double));
    double* dpdY500cyl = (double*) malloc(n * sizeof(double));
    double* dpdY500cyl_eqspaced = (double*) malloc(n * sizeof(double));
    double* normal_Y500cyl = (double*) malloc(n * sizeof(double));
    double* Y500cyl_eqspaced = (double*) malloc(n * sizeof(double));
    double* Y500cyl_distrib = (double*) malloc(n * sizeof(double));
    double* integrand = (double*) malloc(n * sizeof(double));
    double** array_2d = (double**) malloc(n * sizeof(double*));
    Function_1D* fout = (Function_1D*) malloc(sizeof(Function_1D));
    for(i = 0; i < n; i++) array_2d[i] = (double*) malloc(n * sizeof(double));
    fout->x = (double*) malloc(n * sizeof(double));
    fout->y = (double*) malloc(n * sizeof(double));
    
    // Y500cyl(M)
    for (i = 0 ; i < n; i++)
      Y500cyl_vect[i] = m500c2Y500_cyl(p_m_lambda->x[i]*1.e14, redshift, h0, Om, Ol, w0, wa);
    
    // Calculate dpdY500cyl from dpdm
    for(i = 0 ; i < n-1; i++){
      if (i < n-1) dY500cyldm = (Y500cyl_vect[i+1]-Y500cyl_vect[i])/(p_m_lambda->x[i+1]-p_m_lambda->x[i]);
      dpdY500cyl[i] = p_m_lambda->y[i] / dY500cyldm;
    }
    dpdY500cyl[n-1] = 0.;

    // equally-spaced Y500cyl grid AS**** Check range
    double Y500cyl_max = Y500cyl_vect[n-1]*10;
    double Y500cyl_min = Y500cyl_vect[0]/10;
    dY500cyl = (Y500cyl_max - Y500cyl_min) / ((double) n);

    for(i = 0; i < n; i++)                                                                                                     
      {             
	Y500cyl_eqspaced[i] = pow(10,log10(Y500cyl_min)+i*(log10(Y500cyl_max)-log10(Y500cyl_min))/(n-1));  
      }
    
    // AS *** Linear steps
/*     for(i = 0 ; i < n; i++) */
/*       Y500cyl_eqspaced[i] = Y500cyl_min+((double)i)*dY500cyl; */

/* 	if(Y500cyl>Y500cyl_eqspaced[n-1]){ */
/* 	  printf("WARNING: Y500cyl %lf > Y500cyl_eqspaced %lf\n",Y500cyl,Y500cyl_eqspaced[n-1],p_m_lambda->x[n-1]); */
/* 	  printf("WARNING: Y500cyl min %lf, Y500cyl_max %lf\n",Y500cyl_min,Y500cyl_max); */
/* 	  //	  for(i=0;i<n;i++)printf("%lf %lf %lf\n",Y500cyl_eqspaced[i],Y500cyl_vect[i],p_m_lambda->x[i]); */
/* 	} */
/*     if(Y500cyl<Y500cyl_eqspaced[0]) */
/*       printf("WARNING: Y500cyl %lf < Y500cyl_eqspaced %lf\n",Y500cyl,Y500cyl_eqspaced[0],p_m_lambda->x[0]); */
    
	// Interpolate dpdY500cyl to match the values in Y500cyl_eqspaced and put it into dpdY500cyl_eqspaced
    Interp1D_Vector(Y500cyl_vect, dpdY500cyl, n, Y500cyl_eqspaced, dpdY500cyl_eqspaced, n);
/*     for(i=0;i<n;i++){ */
/*       printf("%f %f %f %f \n",Y500cyl_vect[i],dpdY500cyl[i],Y500cyl_eqspaced[i],dpdY500cyl_eqspaced[i]); */
/*     } */


    //AS *** set to zero for neg Y500cyl or to 1? Shouldn't matter?
    for(i=0;i<n;i++)if(Y500cyl_eqspaced[i]<0)dpdY500cyl_eqspaced[i]=0.;
    
    /* Account for the SZ scatter convolve a lognormal distribution with scatter DY500 */
    double Y500cyl_0;
    for (i = 0 ; i < n; i++)
      {
	Y500cyl_0 = Y500cyl_eqspaced[i];
	
        for (j = 0 ; j < n; j++)
	  {
	    lognormal_Y500cyl[j] = 0;
	    if(Y500cyl_0>0)
	      if(Y500cyl_eqspaced[j]>0)
		lognormal_Y500cyl[j] = lognormal(Y500cyl_eqspaced[j], Y500cyl_0, DY500);
	  }
	
        // fix the zero problem
        maxvalue = 1.e-10;
        for(j = 0; j < n; j++){
	  if(lognormal_Y500cyl[j] < maxvalue) lognormal_Y500cyl[j] = maxvalue;
	}
	
        // get the Y500cyl distribution for each initial mass
	//	if(i==n/2)for(j=0;j<n;j++)printf("%d %f %f %f \n",i,Y500cyl_eqspaced[j],normal_Y500cyl[j], lognormal_Y500cyl[j]);		    
	for (j = 0 ; j < n; j++){
	  array_2d[i][j] = lognormal_Y500cyl[j];
	  //	  printf("%d %d %f\n",i,j,array_2d[i][j]);
	}
      }
    
    // Calculate the Y500cyl distribution
    for (i = 0 ; i < n; i++)
      {
        for (j = 0 ; j < n; j++)
	  {
	    integrand[j] = array_2d[j][i] * dpdY500cyl_eqspaced[j];	    
	  }
	//	printf("%f %f \n",integrand[i],dpdY500cyl_eqspaced[i]);	
        Y500cyl_distrib[i] = IntTrap_varstep(Y500cyl_eqspaced, integrand, n, n);
      }
    
    fout->N = n;
    for (i = 0 ; i < n; i++) {
        fout->x[i] = Y500cyl_eqspaced[i];
        fout->y[i] = Y500cyl_distrib[i];
	//	printf("%f %f %f %f \n",p_m_lambda->x[i],p_m_lambda->y[i],Y500cyl_eqspaced[i],Y500cyl_distrib[i]);
    }
    
    free(integrand);
    free(lognormal_Y500cyl);
    free(dpdY500cyl);
    free(dpdY500cyl_eqspaced);
    free(normal_Y500cyl);
    free(Y500cyl_eqspaced);
    for(i = 0; i < n; i++) free(array_2d[i]);
    free(array_2d);
    free(Y500cyl_distrib);    
    return fout;
}
