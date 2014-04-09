/**
 * \file Function.h
 * \brief Function module header
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to handle 1D and 2D functions
 */

#ifndef FUNCTION_2D_H_
#define FUNCTION_2D_H_

#define REGUL_LIN 0
#define REGUL_LOG 1

/**
 * \struct Function_2D 
 * \brief Structure for a 2D function
 */
typedef struct{
  int Nx; /**< int, number of elements in the x (1st) axis */
  int Ny; /**< int, number of elements in the y (2nd) axis */
  double* x; /**< pointer to double, the x vector */
  double* y; /**< pointer to double, the y vector */
  double** z; /**< double pointer to double, the 2D array z(x,y) */
} Function_2D;

/**
 * \struct Function_1D 
 * \brief Structure for a 1D function
 */
typedef struct{
  int N; /**< int, number of elements in the x axis */
  double* x; /**< pointer to double, the x vector */
  double* y; /**< pointer to double, the y(x) vector */
} Function_1D;


Function_2D* Read_Function_2D(char*, char*, char*);

void Copy_Function_1D(Function_1D*, Function_1D**);

short Regularize_1D(Function_1D*, double, short);
void Normalize_1D(Function_1D*);

void Write_2D(Function_2D*, char*);
void Write_1D(Function_1D*, char*);

void Free_Function_2D(Function_2D*);
void Free_Function_1D(Function_1D*);


#endif /* FUNCTION_2D_H_ */
