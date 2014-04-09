#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>


#include "Statistics.h"
#include "Parameters.h"
#include "utils.h"
#include "Function.h"
#include "Data.h"
#include "MassConversion.h"
#include "MassFunction.h"

erf_struct* errfunc;
Data_structure* data;
Function_2D* m200c_2_m500c;
char* home;
char* pwd;




/* get some useful variables */
pwd = getenv("PWD");
home = getenv("HOME");



/* Allocate the tabulated error function */
/* Need it because a priori there is no INVERSE error function in the STL
 * and it is too complicated to calculate it each time */
strcpy(data_path, home);
strcat(data_path, "/dev/ScalingRelationFit/data/erf.dat");
errfunc = InitErf(data_path);


/* Read the data */
strcpy(data_path, pwd);
strcat(data_path, "/dataX.dat");
data = Read_Data(data_path, &Ncl, minz, maxz, 0);
printf("%d clusters\n", Ncl);

/* Read mass conversion factors */
strcpy(data_path, home);
strcat(data_path, "/dev/MassConversion/data");
strcpy(filename_z, data_path);
strcpy(filename_m, data_path);
strcpy(filename, data_path);
strcat(filename_z, "/m200c_to_m500c_z.dat");
strcat(filename_m, "/m200c_to_m500c_m.dat");
strcat(filename, "/m200c_to_m500c_values.dat");
m200c_2_m500c = Read_Function_2D(filename_z, filename_m, filename);
for(i = 0; i < m200c_2_m500c->Ny; i++) m200c_2_m500c->y[i] *= 1.e-14;
printf("%dx%d\tconversion factors (m200c-m500c)\n", m200c_2_m500c->Nx, m200c_2_m500c->Ny);


free(home);
free(pwd);
Free_Function_2D(m200c_2_m500c);
FreeErf(&errfunc);

