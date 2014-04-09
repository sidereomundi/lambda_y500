#include <omp.h>
#include "utils.h"
#include "Data.h"
#include "Likelihood.h"
#include "ScalingRelation.h"
#include "Statistics.h"
#include "linalg.h"
#include "Random.h"
#include "Cosmology.h"
#include "Integration.h"
#include "Interpolation.h"
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

// For each cluster, calculate P(y500cyl|lambda) and return the likelihood
double GetLikelihood_Y500cyl(Parameters_struct* params,
			     Function_1D** vect_p_m_lambda_debiased, 
			     int Ncluster, 
			     Function_1D** vect_p_Y500obs_lambda,
			     Data_structure* data ) 
{
  int i,j,n = vect_p_m_lambda_debiased[0]->N;
 
  //  double probability[Ncluster];
  double* probability = (double*) malloc(Ncluster * sizeof(double)); 
  double loglikelihood = 0.,max_mass=0.;
  Function_1D** p_Y500cyl_lambda;
  double integrand[n];
  double** p_Y500_obs_interp_at_predict = (double*) malloc(Ncluster*sizeof(double));
  double** intg = (double*) malloc(Ncluster*sizeof(double));                                                      
  p_Y500cyl_lambda = (Function_1D**) malloc(Ncluster * sizeof(Function_1D*));
  
  #pragma omp parallel for
  for(i=0;i<Ncluster;i++)
    {
      p_Y500_obs_interp_at_predict[i] = (double*) malloc(n*sizeof(double));
      intg[i] = (double*) malloc(n*sizeof(double));
    }

  #pragma omp parallel for
  for(i = 0; i < Ncluster; i++)
    {
      p_Y500cyl_lambda[i] = Get_P_Y500cyl_lambda(params, vect_p_m_lambda_debiased[i], data[i].redshift);
      Normalize_1D(p_Y500cyl_lambda[i]);	            

      Interp1D_Vector(vect_p_Y500obs_lambda[i]->x,vect_p_Y500obs_lambda[i]->y,vect_p_Y500obs_lambda[i]->N,  
		      p_Y500cyl_lambda[i]->x,p_Y500_obs_interp_at_predict[i],p_Y500cyl_lambda[i]->N); 
      
      for(j=0;j<n;j++)
	{
	  intg[i][j] = p_Y500cyl_lambda[i]->y[j]*p_Y500_obs_interp_at_predict[i][j];
/* 	  printf("%d %f %f %f %f %f %f %f %f \n",i,vect_p_m_lambda_debiased[i]->x[j],vect_p_m_lambda_debiased[i]->y[j], */
/* 		 vect_p_Y500obs_lambda[i]->x[j],vect_p_Y500obs_lambda[i]->y[j], */
/* 		 p_Y500cyl_lambda[i]->x[j],p_Y500cyl_lambda[i]->y[j],p_Y500_obs_interp_at_predict[i][j],intg[i][j]); */
	}
      
      // Get the probability product of the two obs and predicted distributinons

      
      probability[i]=IntTrap_varstep(p_Y500cyl_lambda[i]->x,intg[i], n, n);
      //      printf("Probability for %d: %f \n",i,probability[i]);
      //      printf("%d %f \n",i,probability[i]);
     
      if((probability[i]<=0.) || (probability[i]!=probability[i]))
	{
	  //	  printf("WARNING: probability for %d, params %lf %f %f %f, Y500: %f probability: %f\n",i, params->theta[I_lambda_A], params->theta[I_lambda_B], params->theta[I_lambda_C], params->theta[I_lambda_D], data[i].y0,probability[i]);
	  //				printf("WARNING: probability for %d, Y500 %lf, %e\n",data[i].clus_id, data[i].Y500, probability[i]);
	  probability[i] = 1.e-10;
	}
      
      //printf("%s, likelihood %e\n",data[i].name,probability[i]);
    }  


#pragma omp parallel for
  for(i = 0; i < Ncluster; i++)
    {
      free(p_Y500_obs_interp_at_predict[i]);
      free(intg[i]);
      Free_Function_1D(p_Y500cyl_lambda[i]);
    }

  free(p_Y500_obs_interp_at_predict);
  free(intg);
  free(p_Y500cyl_lambda);

  for(i=0;i<Ncluster;i++)
    loglikelihood += log(probability[i]);
  free(probability);
  return loglikelihood;
}

/* This is the main program
	In the following, the underlying value is lambda_0, the measured value is lambda_obs
	* P(lambda_obs|lambda_0)
	* P(lambda_0|lambda_obs) (not yet debiased)
	* P(M|lambda_obs)_debiased
	* P(Y500cyl|M) */
double Likelihood_Lambdacalib(Parameters_struct* params, 
			      int Ncluster, 
			      Data_structure* data,
			      Function_2D* mass_func,
			      y0_to_Y500cyl_structure* y0_to_Y500cyl){
  
  int i,j;
  FILE *fp;       
  
  double loglikelihood = 0.;

  // Variables for the different distributions
  Function_1D** vect_p_lambda0_lambdaobs;
  Function_1D** vect_p_m_lambda_debiased;
  Function_1D** vect_p_Y500obs_lambda;
  vect_p_lambda0_lambdaobs = (Function_1D**) malloc(Ncluster * sizeof(Function_1D*));
  vect_p_m_lambda_debiased = (Function_1D**) malloc(Ncluster * sizeof(Function_1D*));  
  vect_p_Y500obs_lambda = (Function_1D**) malloc(Ncluster * sizeof(Function_1D*));  
	
  double** plan_lambda_m,plan_Y500cyl_obs_m;

  int Nfilt = params->theta[I_Nfilt];
    
  // lambda scaling parameters
  double Alambda = params->theta[I_lambda_A];
  double Blambda = params->theta[I_lambda_B];
  double Clambda = params->theta[I_lambda_C];
  double Dlambda = params->theta[I_lambda_D];
  
  // Cosmology parameters
  
  double h0 = params->theta[I_H0];
  double Om = params->theta[I_OM];
  double Ol = params->theta[I_OL];
  double w0 = params->theta[I_W0];
  double wa = params->theta[I_WA];
 

  // lambda grid in lambda_vect
  double lambda_max = 300;
  double lambda_min = 1;
  int Nlambda = 50;
  double lambda0 = 1.;
  double dlambda = 10.;
  double lambda_vect[Nlambda];
  double R500c_vect[Nlambda];
  double Y500cyl_max = 0.05;
  double Y500cyl_min = -0.005;

  // Linear steps in Lambda
  /* 	for(i = 0; i < Nlambda; i++) */
  /* 	  lambda_vect[i] = lambda0 + dlambda*i; */
  // Switched to log steps in lambda AS
#pragma omp parallel for
  for(i = 0; i < Nlambda; i++)
    {
      //      lambda_vect[i] = lambda_min+i*(lambda_max-lambda_min)/(Nlambda-1);
      lambda_vect[i] = pow(10,log10(lambda_min)+i*(log10(lambda_max)-log10(lambda_min))/(Nlambda-1));
    }
  // Calculate 2D array with lambda-M relation containing Dlambda scatter and uncertainty

  plan_lambda_m = Get_array_lambda_M(lambda_vect, Nlambda, Dlambda);
  // Calculate 2D array with plan_Y500cyl_obs-M relation containing measurment uncertainty

  // Get P(lambda_0|lambda_obs) (not yet debiased)
#pragma omp parallel for
  for(i = 0; i < Ncluster; i++)
    {
      vect_p_lambda0_lambdaobs[i] = Get_P_lambda0_lambdaobs(plan_lambda_m, lambda_vect, Nlambda, data[i].lambda,data[i].dlambda);
    }

  
#pragma omp parallel for
  for(i = 0; i < Nlambda; i++)
    free(plan_lambda_m[i]);
  free(plan_lambda_m);
  
  /* Get P(M|lambda)_debiased: P(lambda_0|lambda_obs) * mass function
     then check if the distribution makes sense
     finally normalize it */
#pragma omp parallel for
    for(i = 0; i < Ncluster; i++)
    {	      	  
      vect_p_m_lambda_debiased[i] = Get_P_M_lambda(params, vect_p_lambda0_lambdaobs[i], mass_func, data[i].redshift);
      Normalize_1D(vect_p_m_lambda_debiased[i]);      
/*       for (j=0;j<Nlambda;j++){printf(" %f %f %f %f \n",vect_p_lambda0_lambdaobs[i]->x[j],vect_p_lambda0_lambdaobs[i]->y[j], */
/* 				     vect_p_m_lambda_debiased[i]->x[j],vect_p_m_lambda_debiased[i]->y[j]);} */
    }
  
#pragma omp parallel for
  for(i = 0; i < Ncluster; i++)
    {	      	  
      vect_p_Y500obs_lambda[i] = Get_P_Y500cyl_obs_lambda(params,vect_p_m_lambda_debiased[i],data[i].redshift,data[i].y0,data[i].dy0,data[i].filter,Nlambda,y0_to_Y500cyl);
      Normalize_1D(vect_p_Y500obs_lambda[i]);      
      //      for (j=0;j<Nlambda;j++){printf("%e %e  \n",vect_p_Y500obs_lambda[i]->x[j],vect_p_Y500obs_lambda[i]->y[j]);}
    }

  //  return 0;
  // ln-Likelihood of the mass calibration
  loglikelihood = GetLikelihood_Y500cyl(params, vect_p_m_lambda_debiased, Ncluster, vect_p_Y500obs_lambda,data);
  
  // Free everything
#pragma omp parallel for
  for(i = 0; i < Ncluster; i++)
    {
      Free_Function_1D(vect_p_lambda0_lambdaobs[i]);
      Free_Function_1D(vect_p_m_lambda_debiased[i]);
      Free_Function_1D(vect_p_Y500obs_lambda[i]);
    }
  free(vect_p_lambda0_lambdaobs);
  free(vect_p_m_lambda_debiased);
  free(vect_p_Y500obs_lambda);
  return loglikelihood;
}


int call_likelihood_library (int      argc, void *   argv[])
{
  Parameters_struct* params; 
  int Ncluster,Nparams,Nfilt;
  Data_structure* data;
  Function_2D* mass_func;
  int i,j;
  double * loglikelihood;
  y0_to_Y500cyl_structure* y0_to_Y500cyl; 

  Ncluster=*(int *) argv[0];
  Nparams=*(int *) argv[1];
  
  struct  Function_Param{
    int N;
    double THETA[Nparams];
  } *PF_IDL;
  
  PF_IDL=((struct Function_Param *)argv[2]);
  params = malloc(sizeof(Parameters_struct));
  /* number of parameters in the likelihood */
  params->N = PF_IDL[0].N;
  params->theta = malloc(params->N*sizeof(double));
  //#pragma omp parallel for
  for (i=0;i<PF_IDL[0].N;i++)
    params->theta[i] = PF_IDL[0].THETA[i];

  Nfilt = params->theta[I_Nfilt];

  struct  Function_MF{
    int Nx;
    int Ny;
    double X[500];
    double Y[51];
    double Z[51][500];
  } *MF_IDL;

  MF_IDL=((struct Function_MF *)argv[3]);
  mass_func = malloc(sizeof(Function_2D));
  mass_func->Nx = MF_IDL[0].Nx; /* 1st dimension is mass */
  mass_func->Ny = MF_IDL[0].Ny; /* 2nd dimension is redshift */
  mass_func->x = malloc(mass_func->Nx*sizeof(double));
  mass_func->y = malloc(mass_func->Ny*sizeof(double));
  mass_func->z = malloc(mass_func->Nx*sizeof(double*));
  for(i = 0; i < mass_func->Nx; ++i) 
    mass_func->z[i] = malloc(mass_func->Ny*sizeof(double));

  for (i=0;i<MF_IDL[0].Nx;i++)
    mass_func->x[i] = MF_IDL[0].X[i];

  for (i=0;i<MF_IDL[0].Ny;i++)
    mass_func->y[i] = MF_IDL[0].Y[i];
  
  for(i = 0; i < mass_func->Nx; i++)
    for(j = 0; j < mass_func->Ny; j++)
      { mass_func->z[i][j] = MF_IDL[0].Z[j][i];
	//	printf("%d %d %lf \n",i,j,MF_IDL[0].Z[j][i]*1.e10);
      }

  struct  Function_Data{
    int     CLUS_ID;    
    int     FIELD_ID;       
    double   REDSHIFT;       
    double   LAMBDA;         
    double   dLAMBDA;      
    double   y0[Nfilt];             
    double   dy0[Nfilt];          
    double   filter[Nfilt];          
  } *DF_IDL;
  

  DF_IDL=((struct Function_Data *)argv[4]);
  //  for (i=0;i<Nfilt;i++){printf("y0: %e dy0: %e filter %f \n",DF_IDL[0].y0[i],DF_IDL[0].dy0[i],DF_IDL[0].filter[i]);};    
  //  for (id=0;id<Ncluster;id++){printf("%d: lambda %f\n",id,DF_IDL[id].LAMBDA);};
  data = malloc(Ncluster * sizeof (Data_structure));
  
  //#pragma omp parallel for
  for (i=0;i<Ncluster;i++)
    {
      data[i].clus_id    =    DF_IDL[i].CLUS_ID;
      data[i].field_id   =    DF_IDL[i].FIELD_ID;
      data[i].redshift   =    DF_IDL[i].REDSHIFT;
      data[i].lambda	 =    DF_IDL[i].LAMBDA;	 
      data[i].dlambda    =    DF_IDL[i].dLAMBDA;
      for(j=0;j<Nfilt;j++)
	{
	  data[i].y0[j]     =    DF_IDL[i].y0[j];
	  data[i].dy0[j]    =    DF_IDL[i].dy0[j];
	  data[i].filter[j] =    DF_IDL[i].filter[j];
	}
    }
  
  y0_to_Y500cyl = ((struct y0_to_Y500cyl_structure *)argv[5]);

  loglikelihood = (double *) argv[6];
  
  /*   for (i=0;i<10;i++){printf("lambda: %f\n",data[i].lambda);};     */
  //  printf("clus_id %d \n",data[0].clus_id);
  //  for (i=0;i<Nfilt;i++){printf("y0: %e dy0: %e filter %f \n",data[0].y0[i],data[0].dy0[i],data[0].filter[i]);};    
  //  for (i=0;i<Ny500interp;i++){printf("r500: %f y0y500: %f\n",y0_to_Y500cyl->r500[i],y0_to_Y500cyl->y0_to_Y500cyl[i]);};    

  *loglikelihood = Likelihood_Lambdacalib(params,Ncluster,data,mass_func,y0_to_Y500cyl);
  
  free(params->theta);  
  free(params);        
  free(mass_func->x);                                                                                                        
  free(mass_func->y);                                                                                                        
  for(i = 0; i < mass_func->Nx; ++i)                                                                                         
    free(mass_func->z[i]);                                                                                             
  free(mass_func->z);                                                                                                        
  free(mass_func);  
  free(data);

  return 0;
}
