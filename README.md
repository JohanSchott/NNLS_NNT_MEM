# Least Square, Tikhonov, Maximum Entropy Method  (LSTMEM)
This is a program for doing analytical continuations using one of the methods:
- Least Square
- Tikhonov
- Maximum Entropy Metod

The input has to be Matsubara frequency data. Three columns are expected in the input file. 
The first should be the Matsubara frequencies, the second the real part of the function and the third column the imaginary part of the function. 

An control file called `ac.inp` is expected to exist in the current folder. 

The spectral function will be calculated. 

A smeared spectral function is also printed as output. 
It is calculated by taking the spectral function (determined just above the real axis) and using the Hilbert transform to evaluate the spectral function a distance $\delta$ above the real axis ($\frac{-1}{\pi}\text{Im}[G(\omega+i\delta)]$). 

