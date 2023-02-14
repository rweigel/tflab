# TFLAB

A Matlab package for computing transfer functions in the frequency domain using various methods.

# Dependencies

Requires

1. MATLAB Signal Processing Toolbox
2. MATLAB Statistics Toolbox

Tested on MATLAB 2018b and 2021a.

# Getting Started

See [`tflab_demo.m`](tflab_demo.m) for example usage and [`demos/`](demos) for examples where the various options described below are used.

The main function is `tflab.m`.

Unit test may be run by executing `tflab_test`.

# Options

* Time domain windowing using any function available in MATLAB's Signal Processing Toolbox
* Prewhitening
* Multiple regression methods (Ordinary Least Squares, Total Least Squares, Weighted, Robust)
* Computing $Z$ in time domain windows and averaging Z
* Computing DFTs in time domain windows and computing Z based on combined frequency bines.
* Binning of frequencies using log and linear windows. (See `evalfreq()`.)

The functions called for the above may be user-defined.

In addition, this package supports

* Interpolation of estimated $Z$
* Prediction (using interpolated $Z$)
* Computation and display of all aspects of $Z$ calculation and evaluation
* Intercomparison of $Z$ computed using different options
* Comparison of computed $Z$ with those found in EDI files
