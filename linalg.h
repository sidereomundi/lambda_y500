/**
 * \file linalg.h
 * \brief Linear Algebra module header
 * \date on: 9 September 2011
 * \author Gurvan
 *
 * Module for basic linear algebra calculations
 */

#ifndef LINALG_H_
#define LINALG_H_

/**
 * \struct MyMatrix
 * \brief Matrix structure
 */
typedef struct{
  int N; /**< an int, length N of the matrix NxN*/
  double** m; /**< s double pointer to double, 2D matrix array */
} MyMatrix;

double* MatVectProd(MyMatrix*, double*);
MyMatrix* MatMatProd(MyMatrix*, MyMatrix*);
double VectVectProd(int, double*, double*);
void Diagonalize(MyMatrix*, MyMatrix*, double*);
MyMatrix* InvMat(MyMatrix*, double*);
MyMatrix* TransposeMat(MyMatrix*);
double DetMat(MyMatrix*);

#endif

