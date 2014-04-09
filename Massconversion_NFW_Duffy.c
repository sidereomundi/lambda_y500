#include <stdlib.h>
#include <stdio.h>
#include "math.h"

#define EPS 1.e-5

// Often needed
double lnterm(double c){
	return log(1.+c)-c/(1.+c);
}

// Mass-concentration relation from Duffy et al 2008
double Get_c200_M200(double M200, double z){
	return 6.71*pow(M200/2.e12, -0.091)*pow(1.+z, -0.44);
}

// Expression to be minimized by Get_c_DELTA_c200
double rootdiff_c(double c200, double cfind, double ratio){
	return lnterm(c200)/lnterm(cfind) - ratio*pow(c200/cfind, 3.);
}

// Rootfind c(DELTA) as a function of c200
double Get_c_DELTA_c200(double c200, double DELTA){
	double minc = 0.01;
	double maxc = 10.;
	double ratio = 200./DELTA;
	double guess;
	
	// Safety first
	if(rootdiff_c(c200, minc, ratio)>0.)
		printf("ERROR, minc is too large!\n");
	if(rootdiff_c(c200, maxc, ratio)<0.)
		printf("ERROR, maxc is too small!\n");
	
	double diff=10.;
	while((diff<-EPS)||(diff>EPS))
	{
		guess = (minc+maxc)/2.;
		diff = rootdiff_c(c200, guess, ratio);
		if(diff<0.)
			minc = guess;
		else
			maxc = guess;
	}
	return guess;
}

// Expression to be minimized by Get_M200_M_DELTA
double rootdiff_M(double M200, double M_in, double DELTA, double z){
	double c200 = Get_c200_M200(M200, z);
	double c_in = Get_c_DELTA_c200(c200, DELTA);
	
	return M200/M_in - lnterm(c200)/lnterm(c_in);
}

// Get M200 as a function of M_in
double Get_M200_M_DELTA(double M_in, double DELTA, double z){
	double ratio = 200./DELTA;
	
	double Mmin = M_in/ratio/10.;
	double Mmax = M_in/ratio*10.;
	
	double guess;
	
	// Safety first
    if (rootdiff_M(Mmin,M_in,DELTA,z) > 0.)
	{
        printf("!! Mmin is too large !!\n");
        exit(1);
	}
    if (rootdiff_M(Mmax,M_in,DELTA,z) < 0.)
	{
        printf("!! Mmax is too small !!\n");
        exit(1);
	}
	
	
	double diff = 10.;
	while((diff<-EPS)||(diff>EPS))
	{
		guess = (Mmin+Mmax)/2.;
		diff = rootdiff_M(guess, M_in, DELTA, z);
		if(diff<0.)
			Mmin = guess;
		else
			Mmax = guess;
	}
	return guess;	
}
