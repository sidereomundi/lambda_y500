// Interface to Gurvan's dispersion likelihood code

#include "cosmomodule.h"

#include "SZFitC.h"
#include "Parameters.h"
#include "Function.h"
#include "Statistics.h"
#include "MassConversion.h"
#include "MassFunction.h"
#include "Data.h"
#include "Likelihood.h"



double driver_masscalib(int Ncl,
            CAT* cluster,
			double **n,
			 struct COSMO *cosmo,
			 struct SCALING *scaling,
			 Data_structure* data,
			 Function_2D* m200c_2_m500c, 
			 erf_struct* errfunc, 
			 short sigma_flag, 
			 short yx_flag, 
			 short debug_flag)
{
	double likelihood;
	
  /* m200c_2_m500c is the mass conversion table which must be read before
   * errfunc is the tabulated error function which must be read before
   * sigma_flag is 1 if you want to you dispersion in mass calibration
   * yx_flag is 1 if you want to use Yx as well in mass calibration
   * debug_flag is 1 if you want to show additional messages and write some files
   */

  double dm;
  double dz;
  int i, j;
	

  /* Objects needed for the mass calibration likelihood */
  /* fit parameters */
  Parameters_struct* params;
	/* 2D mass function */
  Function_2D* mass_func;

  /* memory allocation */
  params = malloc(sizeof(Parameters_struct));
  /* number of parameters in the likelihood */
  params->N = 19;
  params->theta = malloc(params->N*sizeof(double));

	
  mass_func = malloc(sizeof(Function_2D));
  /* number of redshifts and
   * number of masses */
  mass_func->Nx = NMBIN; /* 1st dimension is mass */
  mass_func->Ny = NZBIN; /* 2nd dimension is redshift */
  mass_func->x = malloc(mass_func->Nx*sizeof(double));
  mass_func->y = malloc(mass_func->Ny*sizeof(double));
  mass_func->z = malloc(mass_func->Nx*sizeof(double*));
	for(i = 0; i < mass_func->Nx; ++i) 
		mass_func->z[i] = malloc(mass_func->Ny*sizeof(double));

/* Copy the mass function */
	
  dm = (MASSMAX-MASSMIN)/(NMBIN-1.);
  dz = (ZMAX-ZMIN)/(NZBIN-1.);
	
  /* 10^m/h in Jiayi's MF,
   * m/1e14 in Gurvan's MF */
	
	for(i = 0; i < mass_func->Nx; i++)
		mass_func->x[i] = pow(10., MASSMIN+dm*((double) i ))/(1.e14 * cosmo->hubble);
	
  /* in Gurvan's MF, z goes increasing
   * in Jiayi's, the other way around : z=ZMAX-i*dz */
	
	for(j = 0; j < mass_func->Ny; j++)
	  mass_func->y[mass_func->Ny-1-j] = ZMAX-dz*((double) j);
	
	for(i = 0; i < mass_func->Nx; i++)
	  for(j = 0; j < mass_func->Ny; j++)
      /* 1st and 2nd dimensions are reversed between Jiayi's and Gurvan's definition */
			mass_func->z[i][mass_func->Ny-1-j] = n[j][i];
	
  /* Jiayi's MF in dN/dlogM,
   * Gurvan's is dN/dM */
/* dN/dM = dN/dlogM * dlogM/dM 
       *       = dN/dlogM * 1/(LN10*M) */
	for(i = 0; i < mass_func->Nx; i++)
		for(j = 0; j < mass_func->Ny; j++)
			mass_func->z[i][j] /= log(10.)*mass_func->x[i];
	
	/*
	FILE *fp;
    fp=fopen("MF_n.dat","w");
	
	for(j=0;j<NMBIN;j++)
    {
        for(i=0;i<NZBIN;i++)
            fprintf(fp,"%e ",n[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    
	
	fp=fopen("driver_mf_x.dat","w");
	
	for(i = 0; i < mass_func->Nx; i++)
		fprintf(fp,"%lf\n",log10((1.e14*mass_func->x[i]*cosmo->hubble)));
	
	fclose(fp);
	
	
	fp=fopen("driver_mf_y.dat","w");
	
	for(i = 0; i < mass_func->Ny; i++)
		fprintf(fp,"%lf\n",mass_func->y[i]);
	fclose(fp);
	
	
	fp=fopen("driver_mf.dat","w");
	
	for(i = 0; i < mass_func->Nx; i++) {
		for(j = 0; j < mass_func->Ny; j++)
			fprintf(fp,"%e ",mass_func->z[i][j]);
		fprintf(fp,"\n");
	}
	
	fclose(fp);
    exit(0);
	*/

	
  /* Copy the fit parameters */
  params->theta[I_SZ_A] = scaling->A_SZ;
  params->theta[I_SZ_B] = scaling->B_SZ;
  params->theta[I_SZ_C] = scaling->C_SZ;
  params->theta[I_SZ_D] = scaling->D_SZ;

  params->theta[I_X_A] = scaling->A_Xray;
  params->theta[I_X_B] = scaling->B_Xray;
  params->theta[I_X_C] = scaling->C_Xray;
  params->theta[I_X_D] = scaling->D_Xray;
	
  params->theta[I_S_A] = scaling->A_disp;
  params->theta[I_S_B] = scaling->B_disp;
  params->theta[I_S_C] = scaling->C_disp;
  //params->theta[I_S_D] = scaling->D_disp;

  params->theta[I_lambda_A] = scaling->Alambda;
  params->theta[I_lambda_B] = scaling->Blambda; 
  params->theta[I_lambda_C] = scaling->Clambda;  
  params->theta[I_lambda_D] = scaling->Dlambda;  
      
  params->theta[I_S_D_0] = scaling->D_disp_0;
  params->theta[I_S_D_N] = scaling->D_disp_N;

  params->theta[I_H0] = cosmo->hubble;
  params->theta[I_OM] = cosmo->Omega_m;
  params->theta[I_OL] = cosmo->Omega_l;
  params->theta[I_W0] = cosmo->w0;
  params->theta[I_WA] = cosmo->wa;
    
    
    j=0;
    for(i=0;i<Ncl;i++)
    {
        if((cluster[i].logYxerror==-1.)&&(cluster[i].veldisp==0.))
	  continue;
		        
        strcpy(data[j].name, cluster[i].name);
	strcpy(data[j].field, cluster[i].field);
	data[j].zeta = cluster[i].zeta;
	data[j].zetaerr = cluster[i].zetaerr;
	data[j].redshift = cluster[i].z;
	data[j].lambda = cluster[i].lambda;
	data[j].lambdaerr = cluster[i].lalmbdaerr;			
	j++;
    }
	
	//for(i=0;i<j;i++)
	  //      printf("%s\t%s\t%lf\t%lf\t%lf\t%lf\t%d\t%lf\n", data[i].name,data[i].field, data[i].xi,data[i].redshift,data[i].yx,data[i].yxerr,data[i].N,data[i].sigma);
		
    
	likelihood = GetLikelihood(params, 
			     j, 
			     data,
			     m200c_2_m500c, 
			     mass_func, 
			     errfunc,
			     sigma_flag, // use dispersion for mass calibration
			     yx_flag, // if you wanna use Yx as well in the mass calibration
				 0,
			     0,
			     debug_flag, // =1 if needed, in order to show messages and write some files
			     0,
				 0, 0.);
				 
				 //printf("mass calibration likelihood %e, log %e\n",likelihood, log(likelihood));

  /* free memory */
  free(params->theta);
  free(params);


	free(mass_func->x);
	free(mass_func->y);
	for(i = 0; i < mass_func->Nx; ++i)
		free(mass_func->z[i]);
	free(mass_func->z);
	free(mass_func);

	if((!isfinite(likelihood))||(isnan(likelihood)))
		printf("ERROR: mass calibration loglikelihood %e\n",likelihood);

	return likelihood;
}
