#ifndef PARAMETERS_H_
#define PARAMETERS_H_

/**
 * \struct Parameters_struct
 * \brief Fit parameter structure
 */
typedef struct{
  int N; /**< an int, number of fit parameters */
  double* theta; /**< a pointer to double, array of fit parameters */
} Parameters_struct;

#define I_SZ_A 0
#define I_SZ_B 1
#define I_SZ_C 2
#define I_SZ_D 3

#define I_S_A 4
#define I_S_B 5
#define I_S_C 6
#define I_S_D 7

#define I_S_D_0 17
#define I_S_D_N 18

#define I_X_A 8
#define I_X_B 9
#define I_X_C 10
#define I_X_D 11

#define I_H0 12
#define I_OM 13
#define I_OL 14
#define I_W0 15
#define I_WA 16

#define I_lambda_A 17
#define I_lambda_B 18
#define I_lambda_C 19
#define I_lambda_D 20

#define I_Y500cyl_D 21
#define I_Nfilt 22


#define Ny500interp 146


#endif
