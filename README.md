# TFLAB

A MATLAB package for computing transfer functions in the frequency domain using various methods with many options.

# Dependencies

Requires

1. MATLAB Signal Processing Toolbox
2. MATLAB Statistics and Machine Learning Toolbox

Tested on MATLAB 2018b, 2021a, and 2022b.

# Getting Started

See [`tflab_demo.m`](tflab_demo.m) for example usage and [`demos/`](demos) for examples where the various options described below are used.

The main function is `tflab.m`.

Unit tests may be run by executing `tflab_test`.

# Options

The functions called for the following may be user-defined.

* Time domain windowing using any function available in MATLAB's Signal Processing Toolbox (see `help window` for a list)
* Prewhitening using `filter()` or `yulewalker()`
* Multiple regression methods (Ordinary Least Squares, Total Least Squares, Weighted, Robust). Some methods have multiple implementations. For example, MATLAB's `robustfit()` can be used or an implementation based on references cited in the MT literature can be used.
* Computing $Z$ for time windows and averaging Z (e.g., "stacking")
* Computing DFTs in time windows and computing Z based on combined frequency bins.
* Binning of frequencies using log and linear windows. (See `evalfreq()`.)

In addition, this package supports

* Interpolation of estimated $Z$ onto a uniform frequency grid
* Prediction (using interpolated $Z$)
* Computation and display of all aspects of $Z$ calculation and evaluation
* Intercomparison of $Z$ computed using different options
* Comparison of `tflab` $Z$s those found in IDE, EDI, and XML files (using the [zread](https://github.com/rweigel/zread) reader)
