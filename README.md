# TFLAB

A MATLAB package for computing, comparing, and visualizing multi-dimensional transfer functions ($Z$s) in the frequency domain using various methods.

**Dependencies**

No additional MATLAB toolboxes are required.

**Optional Dependencies**

1. The MATLAB Signal Processing Toolbox is needed for time domain taper windows other than `rectangualar` and `parzen`.
2. MATLAB Statistics and Machine Learning Toolbox is needed to use the `regress` and `robustfit` regression functions; `tflab` has its own version of robust regression that is used by default.

# Getting Started

The main function is `tflab.m`.

See [`tflab_demo.m`](tflab_demo.m) for example usage and [`demos/`](demos) for examples where the various options described below are used.

Unit tests can be run by executing `tflab_test`.

**Example Usage**



# Options

The functions called for the following may be user-defined.

* Time domain tapering using any function available in MATLAB's Signal Processing Toolbox (e.g., Hamming, see `help window` for a list)
* Prewhitening using MATLAB's `filter()` or `yulewalker()`
* Binning of frequencies using log and linear window; see `evalfreq()`.
* Multiple regression methods including Ordinary Least Squares, Total (Orthogonal) Least Squares, Weighted Least Squares, Robust). Some methods have multiple implementations. For example robust fitting can be performed using `tflab`\'s implementation or MATLAB's `robustfit()`, if available.
* Computing $Z$ for time windows and then averaging to get a composite $Z$ (e.g., "stacking")
* Computing DFTs in time windows and computing $Z$ based on combined frequency bins

In addition, this package supports

* Interpolation of estimated $Z$ onto a uniform frequency grid
* Prediction using an interpolated $Z$
* Calculation of metrics in the time (correlation, coherence, prediction efficiency, mean squared error) and frequency domains (coherence and signal to prediction error).
* Logging of all aspects of the $Z$ calculation
* Visuaal intercomparison of $Z$ computed using different options
* Use of the LEMI MT Fortran library to compute $Z$ (this library is not distributed with `tflab`)
* Comparison of `tflab` $Z$s with
  * LEMI MT's $Z$s and
  * those computed by other algorithms and stored IDE, EDI, and XML files (using the [zread](https://github.com/rweigel/zread) reader)
