/**
 * \file linalg.c
 * \brief Linear Algebra module
 * \date on: 9 September 2011
 * \author Gurvan
 *
 * Module for basic linear algebra calculations
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>

#include "linalg.h"

/**
 * \brief Matrix * vector
 * \param m a pointer to MyMatrix, input matrix
 * \param v a pointer to double, input vector
 * \return a pointer to double, the result of the product
 *
 * Allocates and calculates the product of a matrix by a vector
 */
double* MatVectProd(MyMatrix* m, double* v){
  int i,j;
  double* out = malloc(m->N*sizeof(double));
  for(i = 0; i < m->N; i++){
    out[i] = 0;
    for(j = 0; j < m->N; j++)
      out[i] += m->m[i][j] * v[j];
  }
  return out;
}

/**
 * \brief Matrix * matrix
 * \param m a pointer to MyMatrix, 1st input matrix
 * \param n a pointer to MyMatrix, 2nd input matrix
 * \return a pointer to MyMatrix, the result of the product
 *
 * Allocates and calculates the product of a matrix by a matrix
 */
MyMatrix* MatMatProd(MyMatrix* m, MyMatrix* n){
  int i, j, k;
  MyMatrix* out;

  /* allocate new matrix*/
  out = malloc(sizeof(MyMatrix));
  out->N = m->N;
  out->m = malloc(out->N*sizeof(double*));
  for(i = 0; i < out->N; i++) out->m[i] = malloc(out->N*sizeof(double));

  /* do the multiplication */
  for(i = 0; i < out->N; i++){
    for(j = 0; j < out->N; j++){
	out->m[i][j] = 0.;
      for(k = 0; k < out->N; k++){
	out->m[i][j] += m->m[i][k] * n->m[k][j];
      }
    }
  }

  return out;
}

/**
 * \brief Vector * vector
 * \param n a int, length of the vectors
 * \param u a pointer to double, 1st input vector
 * \param v a pointer to double, 2nd input vector
 * \return a double, the result of the product
 *
 *  Calculates the product of a vector by a vector
 */
double VectVectProd(int n, double* u, double* v){
  int i;
  double out = 0.;
  for(i = 0; i < n; i++) out += u[i]*v[i];
  return out;
}

/**
 * \brief Diagonalization
 * \param m a pointer to MyMatrix, input matrix
 * \param eigvec a pointer to MyMatrix, output eigenvectors array
 * \param eigval a pointer to double, output eigenvalues vector
 *
 *  Calulates and sort the eigenvectors and eigenvalues of a matrix
 */
void Diagonalize(MyMatrix* m, MyMatrix* eigvec, double* eigval){
  
  int i, j;
  double data[m->N*m->N];

  gsl_matrix_view mat
    = gsl_matrix_view_array (data, m->N, m->N);
     
  gsl_vector *eval = gsl_vector_alloc (m->N);
  gsl_matrix *evec = gsl_matrix_alloc (m->N, m->N);
     
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (m->N);

  //gsl_vector_view evec_i;

  /* fill in data */
  for(i = 0; i < m->N; i++)
    for(j = 0; j < m->N; j++)
      data[i*m->N+j] = m->m[i][j];

  gsl_eigen_symmv (&mat.matrix, eval, evec, w);
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

  /* get eigen values and eigen vectors */
  for(i = 0; i < m->N; i++){
    eigval[i] = gsl_vector_get (eval, i);
    //evec_i = gsl_matrix_column (evec, i);
    //gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
    for(j = 0; j < m->N; j++) eigvec->m[i][j] = evec->data[j+m->N*i];
  }
     
  gsl_vector_free (eval);
  gsl_matrix_free (evec);

}

/**
 * \brief Matrix inversion
 * \param m a pointer to MyMatrix, input matrix
 * \param det a pointer to double, the output determinant
 * \return a pointer to MyMatrix, inverted  matrix
 *
 *  Allocates and calculates the inverse of the input matrix
 */
MyMatrix* InvMat(MyMatrix* m, double* det){

  int i, j;
  double data[m->N*m->N];

  gsl_matrix_view mat;   
  int s;
  gsl_permutation * p;
  gsl_matrix * inverse;

  MyMatrix* out;

  /* fill in the input matrix */
  for(i = 0; i < m->N; i++)
    for(j = 0; j < m->N; j++)
      data[i*m->N+j] = m->m[i][j];

  mat = gsl_matrix_view_array (data, m->N, m->N); 
  p = gsl_permutation_alloc (m->N);
  inverse = gsl_matrix_alloc(m->N, m->N);

  /* LU decomposition and inversion */
  gsl_linalg_LU_decomp (&mat.matrix, p, &s);
  gsl_linalg_LU_invert (&mat.matrix, p, inverse);
  *det = gsl_linalg_LU_det (&mat.matrix, s);

  /* fill in output matrix */
  out = malloc(sizeof(MyMatrix));
  out->m = malloc(m->N*sizeof(double*));
  out->N = m->N;
  for(i = 0; i < m->N; i++) out->m[i] = malloc(out->N*sizeof(double));

  for(i = 0; i < m->N; i++)
    for(j = 0; j < m->N; j++)
      out->m[i][j] = inverse->data[i*m->N+j];

  gsl_permutation_free (p);
  gsl_matrix_free(inverse);

  return out;

}

/**
 * \brief Matrix determinant
 * \param m a pointer to MyMatrix, input matrix
 * \return a double, determinant
 *
 *  Calculates the determinant of the input matrix
 */
double DetMat(MyMatrix* m){

  int i, j;
  double data[m->N*m->N];

  gsl_matrix_view mat;
  int s;
  gsl_permutation * p;

  double out;

  /* fill in the input matrix */
  for(i = 0; i < m->N; i++)
    for(j = 0; j < m->N; j++)
      data[i*m->N+j] = m->m[i][j];

  /* LU decomposition */
  mat = gsl_matrix_view_array (data, m->N, m->N); 
  p = gsl_permutation_alloc (m->N);

  gsl_linalg_LU_decomp (&mat.matrix, p, &s);
  out = gsl_linalg_LU_det (&mat.matrix, s);

  gsl_permutation_free (p);

  return out;

}

/**
 * \brief Matrix transposition
 * \param m a pointer to MyMatrix, input matrix
 * \return a pointer to MyMatrix, transposed matrix
 *
 *  Allocates and calculates the transposition of the input matrix
 */
MyMatrix* TransposeMat(MyMatrix* m){
  int i, j;
  MyMatrix* out;

  /* allocate new matrix*/
  out = malloc(sizeof(MyMatrix));
  out->N = m->N;
  
  out->m = malloc(out->N*sizeof(double*));
  for(i = 0; i < out->N; i++) out->m[i] = malloc(out->N*sizeof(double));
  /* transpose */
  for(i = 0; i < out->N; i++){
    for(j = 0; j < out->N; j++){
      out->m[i][j] = m->m[j][i];
    }
  }

  return out;
}


