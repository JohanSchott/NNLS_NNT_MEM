# NNLS, NNT, MEM
This is a program for doing analytical continuations using one of the methods:
- Non Negative Least Squares
- Non Negative Tikhonov
- Maximum Entropy Metod

The input has to be Matsubara frequency data. Three columns are expected in the input file. 
The first should be the Matsubara frequencies, the second the real part of the function and the third column the imaginary part of the function. 

A control file called `ac.inp` is expected to exist in the current folder. 

The spectral function will be calculated. 
A smeared spectral function is also printed as output. 
It is calculated by taking the spectral function (determined just above the real axis) and using the Hilbert transform to evaluate the spectral function a distance $\delta$ above the real axis ($\frac{-1}{\pi}\text{Im}[G(\omega+i\delta)]$). 

Compile by copying the example Makefile: `src/Makefile_example` to `src/Makefile` and adjust it to fit to the current machine. In the `src` folder, run: `make` and the binary `NNLS_NNT_MEM` should be created.
Run the program by executing the binary.
