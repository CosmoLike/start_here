all:
	gcc -O3 -shared -o like_fourier.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -funroll-loops -std=gnu99

