CC=gcc 

CFLAGS = -g -fopenmp -m64 -O3  -Wall -g3 -fPIC  -ld -shared  $(GSL_INC) -I /home/moon/saro/LAMBDA_MASS_CALIBRATION/SB/trunk/SebCosmo
LDFLAGS= -fopenmp -Wall -g3 -lm -lgsl -lgslcblas -L /home/moon/saro/LAMBDA_MASS_CALIBRATION/SB/trunk/SebCosmo -lSebCosmo -lgomp


OBJECTS=Convolution.o Cosmology.o Function.o Integration.o Interpolation.o Likelihood.o MassConversion.o MassFunction.o Parameters.o PMC.o Random.o ScalingRelation.o Statistics.o SZFitC.o utils.o linalg.o MixtureModel.o linalg.o

lib: 	Convolution.o Function.o Integration.o Interpolation.o Likelihood.o MassConversion.o Random.o ScalingRelation.o Statistics.o utils.o linalg.o MixtureModel.o 
	$(CC) -shared $(LDFLAGS) -Wl,-soname,Lambdacalib.so -o Lambdacalib.so $^

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $< 

clean:
	rm *.o *.so
