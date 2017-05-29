#!/bin/bash 

rm ./like_fourier

gcc -Wno-missing-braces -Wno-missing-field-initializers -I/usr/local/include -L/usr/local/lib  -o ./like_fourier like_fourier.c -lfftw3 -lgsl -lgslcblas -lm 

./like_fourier

