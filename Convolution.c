/**
 * \file Convolution.c
 * \brief Convolution module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module to calculate the convolution of a function by a kernel
 */

#include <stdio.h>

#include "utils.h"
#include "Convolution.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

/**
 * \brief Makes the convolution f*g = fg
 * \param f a pointer to double
 * \param g a pointer to double
 * \param len_f an int, the length of the array f
 * \param len_g an int, the length of the array g
 * \param fg a pointer to double, result of the convolution
 *
 * Checks that the lengths of f and g are equal and returns the convolution of f by g on the same grid
 */
void Convolve(double* f, double* g, int len_f, int len_g, double* fg){

  int i;
  int n = len_f;
  int cpt;
  
  if (n != len_g) {
    printf("Convolution::Convolve: size(f) != size(g), return f\n");
    for (i = 0; i < n; i++) fg[i] = f[i];
  }
  
  // size *2
  double f_d[2*n];
  double g_d[2*n];
  // initialization
  for (i = 0; i < 2*n ; i++){
    f_d[i] = 1.e-6;
    g_d[i] = 1.e-6;
  }
  for (i = 0; i < n ; i++) f_d[i+n/2] = f[i];
  
  for (i = 0; i < n ; i++) g_d[i+n/2] = g[i];
  
  n *= 2;
  
  // shift g
  double g_s[n];
  cpt = 0;
  for (i = n/2; i < n; i++) g_s[cpt++] = g_d[i];
  for (i = 0; i < n/2 ; i++) g_s[cpt++] = g_d[i];
  
  gsl_fft_complex_wavetable * wave_f;
  gsl_fft_complex_wavetable * wave_g;
  gsl_fft_complex_workspace * work_f;
  gsl_fft_complex_workspace * work_g;
  
  work_f = gsl_fft_complex_workspace_alloc (n);
  work_g = gsl_fft_complex_workspace_alloc (n);
  wave_f = gsl_fft_complex_wavetable_alloc (n);
  wave_g = gsl_fft_complex_wavetable_alloc (n);
  
  double data_f[2*n];
  double data_g[2*n];
  for (i = 0; i < n ; i++){
    data_f[2*i] = f_d[i];
    data_g[2*i] = g_s[i];
    data_f[2*i+1] = 0.;
    data_g[2*i+1] = 0.;
  }
  
  gsl_fft_complex_transform(data_f, 1, n, wave_f, work_f, (gsl_fft_direction)-1);
  gsl_fft_complex_transform(data_g, 1, n, wave_g, work_g, (gsl_fft_direction)-1);
  
  gsl_fft_complex_workspace_free (work_f);
  gsl_fft_complex_workspace_free (work_g);
  gsl_fft_complex_wavetable_free (wave_f);
  gsl_fft_complex_wavetable_free (wave_g);
  
  // Multiply in Fourier space
  double data_fg[2*n];
  for (i = 0; i < n ; i++){
    data_fg[2*i] = data_f[2*i] * data_g[2*i] - data_f[2*i+1] * data_g[2*i+1];
    data_fg[2*i+1] = data_f[2*i] * data_g[2*i+1] + data_f[2*i+1] * data_g[2*i];
  }
  
  gsl_fft_complex_workspace* work_h = gsl_fft_complex_workspace_alloc (n);
  gsl_fft_complex_wavetable* wave_h = gsl_fft_complex_wavetable_alloc (n);
  
  gsl_fft_complex_inverse (data_fg, 1, n, wave_h, work_h);
  
  gsl_fft_complex_workspace_free (work_h);
  gsl_fft_complex_wavetable_free (wave_h);
  
  n /= 2;
  for (i = 0; i < n; i++){
    fg[i] = data_fg[2*(i+n/2)];
  }

}

