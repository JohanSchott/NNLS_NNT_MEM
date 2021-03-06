# NNLS, NNT, MEM
## Goal
This is a program for doing analytical continuations using one of the methods:
- Non Negative Least Squares
- Non Negative Tikhonov
- Maximum Entropy Metod

## Compile
Compile by copying the example Makefile: `src/Makefile_example` to `src/Makefile` and adjust it to fit to the current machine.
In the `src` folder, run: `make` and the binary `NNLS_NNT_MEM` should be created.
LAPACK is required. 
Quadruple precision modified LAPACK routines are stored in the folder `src/quad`.

## How to use it
Run the program by executing the binary.
A control file called `ac.inp` is expected to exist in the current folder. 
Among many things is the Matsubara input filename specified in `ac.inp`.
The input Matsubara data has to be Matsubara frequency data. Three columns are expected in the input file. 
The first should be the Matsubara frequencies, the second the real part of the function and the third column the imaginary part of the function. 
The `ac.inp` is row format strict, one can not change the order of the parameters.
A typical `ac.inp` file exist in folder `tests`. Copy and adjust it for your analtyical continuation. 

### Tests
In the `tests` folder a few test models exist:
- `betheU0`
- `betheU4`
- `Haverkort_wc1_dw0.5`
- `Sm7`

A good start for using the program is to try these tests. Perhaps a good beginning is also to modify the attached `ac.inp` files. One thing one can do is to try both NNLS, NNT and MEM for the test models.

### Output
The spectral function will be calculated. 
A smeared spectral function is also printed as output. 
It is calculated by taking the spectral function (determined just above the real axis) and using the Hilbert transform to evaluate the spectral function a distance $\delta$ above the real axis ( $\frac{-1}{\pi}\text{Im}[G(\omega+i\delta)]$ ). 

## Possible improvements
- Implement kernels for imaginary time
- Test MINPACK's non-linear equation solver
- Implement stochastic sampling methods  
