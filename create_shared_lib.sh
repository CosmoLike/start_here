#!/bin/bash 

rm like_fourier.so

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib -shared -o like_fourier.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm



