/**
 * \file MassConversion.c
 * \brief MassConversion module
 * \date on: 8 August 2011
 * \author Gurvan
 *
 * Module for conversion between mass definitions
 */

#include <stdlib.h>
#include <stdio.h>
#include "math.h"

#include "MassConversion.h"
#include "Interpolation.h"

/**
 * \brief Makes mass conversion
 * \param f a pointer to Function_2D, mass conversion table
 * \param mass_in a double, the input mass
 * \param redshift a double, the input redshift
 * \return a double, the converted mass
 *
 * Makes mass conversion given the input mass, redshift and conversion table
 */
double MassConversion(Function_2D* f, double mass_in, double redshift){
	double factor = Interp2D_Scalar(f->x, f->y, f->z, f->Nx, f->Ny, redshift, mass_in);
	return mass_in * factor;
}

/**
 * \brief Intermediate function
 * \param convert_x a pointer to double, the conversion table masses
 * \param convert_mass a pointer to double, the conversion table interpolated in redshift for each mass
 * \param len a int, length of the conversion table
 * \param mass_in a pointer to double, the input mass vector
 * \param mass_out a pointer to double, the output converted mass vector
 * \param len a int, length of the mass vector
 *
 * Makes the interpolation in mass for the 2D mass conversion (redshift, mass)
 */
void MassConversion_optimal1D(double* convert_x, double* convert_mass, int len, double* mass_in, double* mass_out, int len_mass){
	int i;
	double factor_out[len_mass];
	Interp1D_Vector(convert_x, convert_mass, len, mass_in, factor_out, len_mass);
	for(i = 0; i< len_mass; i++){
		mass_out[i] = mass_in[i] * factor_out[i];
	}
}

/**
 * \brief Makes mass conversion for a vector of input masses at the same redshift
 * \param f a pointer to Function_2D, mass conversion table
 * \param mass_in_vect a pointer to double, the input mass vector
 * \param mass_in_len a int, length of the mass vector
 * \param redshift a double, the single input redshift
 * \param mass_out_vect a pointer to double, the output converted mass vector
 *
 * Makes the optimal mass conversion given the input masses at the samne redshift and conversion table
 */
void MassConversion_optimal2D(Function_2D* f, double* mass_in_vect, int len_mass_in, double redshift, double* mass_out_vect){
	int i, j;
	double convert_tmp[f->Nx];
	double convert_mass[f->Ny];
	for(i = 0; i < f->Ny; i++){
		/* first get the conversion vector for each mass */
		for(j = 0; j < f->Nx; j++) convert_tmp[j] = f->z[j][i];
		/* interpolate at the actual redshift value */
		convert_mass[i] = Interp1D_Scalar(f->x, convert_tmp, f->Nx, redshift);
	}
	MassConversion_optimal1D(f->y, convert_mass, f->Ny, mass_in_vect, mass_out_vect, len_mass_in);
}

/**
 * \brief Makes inverse mass conversion
 * \param f a pointer to Function_2D, mass conversion table
 * \param mass_in a double, the input mass
 * \param redshift a double, the input redshift
 * \return a double, the output converted mass
 *
 * Makes the inverse mass conversion given the input mass, redshift and conversion table. The function is recursive.
 */
double InvMassConversion(Function_2D* f, double mass_in, double redshift, double mass_tmp, double epsilon){
	double factor = Interp2D_Scalar(f->x, f->y, f->z, f->Nx, f->Ny, redshift, mass_tmp);
	double mass_out = mass_in / factor;
	if (fabs(mass_tmp-mass_out) < epsilon)
		return mass_out;
	else
		return InvMassConversion(f, mass_in, redshift, mass_out, epsilon);
}

/**
 * \brief Intermediate function
 * \param convert_x a pointer to double, the conversion table masses
 * \param convert_mass a pointer to double, the conversion table interpolated in redshift for each mass
 * \param len a int, length of the conversion table
 * \param mass_in a pointer to double, the input mass vector
 * \param mass_out a pointer to double, the output converted mass vector
 * \param len a int, length of the mass vector
 *
 * Makes the interpolation in mass for the 2D inverse mass conversion (redshift, mass). The function is recursive.
 */
double InvMassConversion_optimal1D(double* convert_x, double* convert_mass, int len, double mass_in, double mass_tmp, double epsilon){
	double factor = Interp1D_Scalar(convert_x, convert_mass, len, mass_tmp);
	double mass_out = mass_in / factor;
	if (mass_out<convert_x[0]) return mass_out;
	if (mass_out>convert_x[len-1]) return mass_out;
	if (fabs(mass_tmp-mass_out) < epsilon)
		return mass_out;
	else
		return InvMassConversion_optimal1D(convert_x, convert_mass, len, mass_in, mass_out, epsilon);
}

/**
 * \brief Makes inverse mass conversion for a vector of input masses at the same redshift
 * \param f a pointer to Function_2D, mass conversion table
 * \param mass_in_vect a pointer to double, the input mass vector
 * \param mass_in_len a int, length of the mass vector
 * \param redshift a double, the single input redshift
 * \param mass_out_vect a pointer to double, the output converted mass vector
 *
 * Makes the optimal inverse mass conversion given the input masses at the samne redshift and conversion table.
 */
void InvMassConversion_optimal2D(Function_2D* f, double* mass_in_vect, int len_mass_in, double redshift, double epsilon, double* mass_out_vect){
	int i, j;
	double convert_tmp[f->Nx];
	double convert_mass[f->Ny];

	for(i = 0; i < f->Ny; i++){
		/* first get the conversion vector for each mass */
		for(j = 0; j < f->Nx; j++) convert_tmp[j] = f->z[j][i];
		/* interpolate at the actual redshift value */
		convert_mass[i] = Interp1D_Scalar(f->x, convert_tmp, f->Nx, redshift);
	}
	/* loop over input masses */
	for(i = 0; i < len_mass_in; i++){
		mass_out_vect[i] = InvMassConversion_optimal1D(f->y, convert_mass, f->Ny, mass_in_vect[i], mass_in_vect[i]/convert_mass[f->Ny/2], epsilon);
	}
}




