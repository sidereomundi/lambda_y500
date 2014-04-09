/**
 * \file Function.c
 * \brief Function module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to handle 1D and 2D functions
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "Integration.h"
#include "Interpolation.h"
#include "Function.h"

/**
 * \brief Copy a 1D function
 * \param in a pointer to Function_1D, input function
 * \param out a double pointer to Function_1D, output function
 *
 * Allocate a Function_1D structure and make a copy
 */
void Copy_Function_1D(Function_1D* in, Function_1D** out){

  int i;

  /* allocate memory for output */
  *out = (Function_1D*) malloc(sizeof(Function_1D));
  (*out)->N = in->N;
  (*out)->x = (double*) malloc((*out)->N * sizeof(double));
  (*out)->y = (double*) malloc((*out)->N * sizeof(double));

  /* copy */
  for(i = 0; i < (*out)->N; i++){
    (*out)->x[i] = in->x[i];
    (*out)->y[i] = in->y[i];
  }

}

/**
 * \brief Read a 2D function
 * \param filename_m a pointer to char, path of the input file containing the x axis
 * \param filename_z a pointer to char, path of the input file containing the y axis
 * \param filename a pointer to char, path of the input file containing the 2D array
 *
 * Read a 2D function from input files, allocate a Function_2D structure, fill it and return it.
 */
Function_2D* Read_Function_2D(char* filename_m, char* filename_z, char* filename){
  FILE* myfile;
  FILE* myfile_m;
  FILE* myfile_z;
  Function_2D* massconv;
  char buffer[8192];
  int i, j;
  int n;
  int n_m;
  int n_z;
  double value;
  double m_vect[1000];
  double z_vect[1000];
  double** f_matrix = (double**) malloc(1000 * sizeof(double*));
  for(i = 0; i < 1000; i++) f_matrix[i] = (double*) malloc(1000 * sizeof(double));

  /* Read mass axis */
  n_m = 0;
  myfile_m = fopen(filename_m, "r");
  while (fgets(buffer, sizeof(buffer), myfile_m) != 0){
    n = sscanf(buffer, "%lf", &value);
    if (n != 1)
      printf("Failed to read %s\n", filename_m);
    m_vect[n_m++] = value;
  }
  fclose(myfile_m);

  /* Read redshift axis */
  n_z = 0;
  myfile_z = fopen(filename_z, "r");
  while (fgets(buffer, sizeof(buffer), myfile_z) != 0){
    n = sscanf(buffer, "%lf", &value);
    if (n != 1)
      printf("Failed to read %s\n", filename_z);
    z_vect[n_z++] = value;
  }
  fclose(myfile_z);

  /* Read 2D array */
  myfile = fopen(filename, "r");
  //while (fgets(buffer, sizeof(buffer), myfile) != 0){
  for(i = 0; i < n_m; i++){
    for(j = 0; j < n_z; j++){
      fscanf(myfile, "%lf", &f_matrix[i][j]);
      //printf("%lf\t", f_matrix[i][j]);
    }
    //printf("\n");
  }
  fclose(myfile);

  massconv = (Function_2D*) malloc(sizeof(Function_2D));
  /* allocate memory */
  massconv->Nx = n_m;
  massconv->Ny = n_z;
  massconv->x = (double*) malloc(n_m * sizeof (double));
  massconv->y = (double*) malloc(n_z * sizeof (double));
  massconv->z = (double**) malloc(n_m * sizeof (double*));

  for(i = 0; i < n_m; i++){
    massconv->x[i] = m_vect[i];
  }
  for(i = 0; i < n_z; i++){
    massconv->y[i] = z_vect[i];
  }
  for(i = 0; i < n_m; i++){
    massconv->z[i] = (double*) malloc(n_z * sizeof (double));
    for(j = 0; j < n_z; j++){
      massconv->z[i][j] = f_matrix[i][j];
    }
  }

  for(i = 0; i < 1000; i++) free(f_matrix[i]);
  free(f_matrix);

  return massconv;
}

/**
 * \brief Fix some particular problems in a 1D function
 * \param f a pointer to Function_1D, input function to check
 * \param location a double, maximum allowed location of the maximum of the distribution
 * \param regul a short, 0 for linear scale, 1 for log scale, used for symetrization of the function if needed
 * \return a short, 1 if successful, 0 if not
 *
 * Detect and fix problems related to 1D distributions: check for positive values, find a maximum, check that the distribution is a bijection before the location of the maximum, try to correct if not (symmetrizes the distribution in log scale).
 */
short Regularize_1D(Function_1D* f, double location, short regul){
  int i, j;
  int imax;
  double max;
  int imin;
  double min;
  int step = 1;
  double maxfrac = 0.02;
  double minx = f->x[0];
  double maxx = f->x[f->N-1];
  double* interp;

    if(location != 0){
    maxx = location;
  }

  /* find a maximum */
  imax = 0;
  max = -1.e10;
  for(i = step; i < f->N-step; i++){
    if(f->x[i] >minx && f->x[i]<maxx && f->y[i]>max && f->y[i]>f->y[i+step] && f->y[i]>f->y[i-step]){
      max = f->y[i];
      imax = i;
    }
  }

  /* if no maximum in that range, error */
    if(imax == 0) {
        for(i=1;i<f->N;i++)
            if(f->y[0] <= f->y[i]) {
                printf("No maximum found. max %e, i=1 %e\n",max,f->y[1]);
                return 0;
            }
        max = f->y[0];
        //return 0;
    }

  /* find a minimum between 0 and imax */
  imin = 0;
  min = max;
  for(i = 0; i < imax; i++){
    if(f->y[i]<min){
      min = f->y[i];
      imin = i;
    }
  }

  /* if numerical issue */
  if(min < 0. && fabs(min)>.1*max){
    /*FILE* checkfile = fopen("error.dat", "w");
      for(i = 0; i < f->N; i++) fprintf(checkfile, "%lf\t%lf\n", f->x[i], f->y[i]);
      fclose(checkfile);*/
      printf("Numerical issue in regularize. min %lf\n",min);
    return 0;
  }

  /* if too high y for low x values... numerical issue */
  if(min > maxfrac * max && imax!=0){
      
    //      printf("min>maxfrac*max. min %e imin %d max %e imax %d\n",min,imin,max,imax);
    /* just multiply by a polynome of order 1 */
    /*
      for(i = imin; i < imax; i++){
      f->y[i] *= ((float) (i-imin)) / ((float) (imax-imin));
      }
      for(i = 0; i < imin; i++){
      f->y[i] = 0.;
      }
    */
    /* better to symetrize with the right-hand side */
    if(regul == REGUL_LIN){ /* linear scale */
      for(i = 0; i < imax; i++){
          j = 2*imax - i;
          if(j<f->N){
              f->y[i] = f->y[j];
          }else{
              f->y[i] = 0.;
          }
      }
    }
    else{ /* log scale */
      interp = (double*) malloc((imax)*sizeof(double));
      j = 0;
      for(i = 0; i < imax; i++)
          interp[j++] = Interp1D_Scalar(f->x, f->y, f->N, exp( 2.*log(f->x[imax]) - log(f->x[i]) ) );

      j = 0;
      for(i = 0; i < imax; i++)
          f->y[i] = interp[j++];
      
      free(interp);
    }
  }
  else{
    /* Anyway put zero below xmin*/
    for(i = 0; i < imin; i++){
      f->y[i] = 0.;
    }
  }

    /* regularization is ok */
  return 1;

}

/**
 * \brief Normalize a 1D function
 * \param f a pointer to Function_1D, input function to normalize
 *
 * Normalize a 1D function so that the integral is 1
 */
void Normalize_1D(Function_1D* f){
  int i;
  double norm;
  norm = IntTrap_varstep(f->x, f->y, f->N, f->N);
  for(i = 0; i < f->N; i++) f->y[i] /= norm;
}

/**
 * \brief Write a 2D function
 * \param f a pointer to Function_2D, input function to write
 * \param filename a pointer to char, output file path
 *
 * Write a 2D function to a file
 */
void Write_2D(Function_2D* f, char* filename){
  int i, j;
  FILE* myfile = fopen(filename, "w");
  for(i = 0; i < f->Nx; i++){
    for(j = 0; j < f->Ny; j++){
      fprintf(myfile, "%e\t%e\t%e\n", f->x[i], f->y[j], f->z[i][j]);
    }
  }
  fclose(myfile);
}

/**
 * \brief Write a 1D function
 * \param f a pointer to Function_1D, input function to write
 * \param filename a pointer to char, output file path
 *
 * Write a 1D function to a file
 */
void Write_1D(Function_1D* f, char* filename){
  int i;
  FILE* myfile = fopen(filename, "w");
  for(i = 0; i < f->N; i++){
    fprintf(myfile, "%e\t%e\n", f->x[i], f->y[i]);
  }
  fclose(myfile);
}

/**
 * \brief Deallocate a 2D function
 * \param f a pointer to Function_2D, input function to deallocate
 *
 * Deallocate a 2D function
 */
void Free_Function_2D(Function_2D* f){
  int i;
  int n = f->Nx;
  free(f->x);
  free(f->y);
  for(i = 0; i < n; i++) free(f->z[i]);
  free(f->z);
  free(f);
}

/**
 * \brief Deallocate a 1D function
 * \param f a pointer to Function_1D, input function to deallocate
 *
 * Deallocate a 1D function
 */
void Free_Function_1D(Function_1D* f){
  free(f->x);
  free(f->y);
  free(f);
}


