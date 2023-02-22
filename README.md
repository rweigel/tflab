# TFLAB

A MATLAB package for computing, comparing, and visualizing multi-dimensional transfer functions in the frequency domain using various methods.

# Optional Dependencies

1. MATLAB Signal Processing Toolbox (to use time domain windows other than `rectangualar` and `Parzen`.)
2. MATLAB Statistics and Machine Learning Toolbox for robust regression; this package has its own version of robust regression that can be used instead.

# Getting Started

See [`tflab_demo.m`](tflab_demo.m) for example usage and [`demos/`](demos) for examples where the various options described below are used.

The main function is `tflab.m`.

Unit tests may be run by executing `tflab_test`.

# Options

The functions called for the following may be user-defined.

* Time domain windowing using any function available in MATLAB's Signal Processing Toolbox (e.g., Hamming, see `help window` for a list)
* Prewhitening using `filter()` or `yulewalker()`
* Multiple regression methods (Ordinary Least Squares, Total Least Squares, Weighted Least Squares, Robust Regression). Some methods have multiple implementations. For example, MATLAB's `robustfit()` can be used or an implementation based on references cited in the MT literature can be used.
* Computing $Z$ for time windows and then averaging to get composite $Z$ (e.g., "stacking")
* Computing DFTs in time windows and computing $Z$ based on combined frequency bins
* Binning of frequencies using log and linear window; see `evalfreq()`.

In addition, this package supports

* Interpolation of estimated $Z$ onto a uniform frequency grid
* Prediction (using interpolated $Z$)
* Calculation of metrics in the time and frequency domains: Correlation, Coherence, Prediction Efficiency, Mean Squared Error, and Signal to Prediction Error.
* Computation and display of all aspects of $Z$ calculation and evaluation
* Intercomparison of $Z$ computed using different options
* Use of the LEMI MT Fortran library to compute $Z$ (this library is not distributed with `tflab`)
* Comparison of `tflab` $Z$s with LEMI MTs and those found in IDE, EDI, and XML files (using the [zread](https://github.com/rweigel/zread) reader)
