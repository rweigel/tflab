# transferfn
Estimate transfer functions from data

In MT data analysis, one seeks to find tf of form

(1) Ex = ZxxBx + ZxyBy
(2) Ey = ZyxBx + ZyyBy

Given time domain measurements of E and B.

(1) and (2) can be treated separately. The 1-D case of E = ZB has been treated extensively in the literature.

Due to noise, we need to specify how Z was solved for. For example, minimize SSE of E(t)-E_m(t).

There are many complications that one must account for when chosing a method:

1. Efficiency of calculation
2. Non-stationary noise
3. Non-Gaussian noise including spikes and baseline offset changes

Most MT methods for computing Z are done in the frequency domain (using the FFT to transform the data from the time to frequency domain),

In this case, case one must account for
1. How non-stationary noise in time domain manifests in the fd

2. How non-gaussian noise, including spikes and baseline offsets, in time domain manifests in the fd - a single-point spike in the time domain affects all frequencies in the fd.

3. The spectrum of the measurements are not white. Most statistical results are for a flat spectrum. Transforming the spectra to be flat introduces an error that must be understood. See for example Austin 1998 for comparsions of estimating power law slopes w/ and w/o prewhitening.

4. Aliasing

5. Spectral leakage - 

6. Issues associated with regression when the regression error is not Gaussian ("robustness")

7. The location and width of frequency bands used for regression

8. Data gaps

In the MT literature there is a wide variation in the assumptions and methods used.

Existing code is [1] and [2] and both have been extensively used in the literature. [1] Is available on-line for use. [2] is available by request. Code [3] is shipped with MT instruments.