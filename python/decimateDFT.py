import numpy as np
from numpy.fft import fft, fftfreq
from scipy import signal

def decimateDFT(x, n, factor=2, all=False):

    # TODO: n and factor must be int.
    # len(x) >= 27 (for signal.decimate, which requires padlen >= 27)

    if all:
        k = 0
        X = []
        while x.shape[0] >= n*factor:
            (fo, x_dft_aveo, x) = decimateDFT(x, n, factor=2, all=False)
            if k == 0:
                f = fo
                x_dft_ave = x_dft_aveo
            else:
                f = np.vstack((f, fo/np.power(factor, k)))
                x_dft_ave = np.vstack((x_dft_ave, x_dft_aveo))
            X.append(x)
            k = k + 1

        return (f, x_dft_ave, X)
    
    if False:
        x = signal.decimate(x, factor)

    if True:
        from numpy.fft import ifft
        xf = fft(x)
        fx = fftfreq(x.shape[0])
        I = np.abs(fx) > 0.25
        xf[I] = 0
        #print(np.std(x))
        so = np.std(x)
        print(x.shape[0])
        print('std before lpf = %.3f' % so)
        x = ifft(xf)
        sf = np.std(x)
        print('std after lpf = %.3f' % sf)
        print('Ratio %.2f' % (so/sf))
        x = x[::factor]
        print('std after decimation = %.3f' % np.std(x))
        
    if False:
        # Code of signal.decimate used in above
        system = signal.dlti(*signal.cheby1(8, 0.05, 0.8 / factor))
        b, a = system.num, system.den
        x = signal.filtfilt(b, a, x) # Zero phase filter
        #x = signal.lfilter(b, a, x)
        x = x[::2]
        #w, h = signal.freqz(b, a)
        #plt.semilogy(w/np.pi, abs(h), 'b')
    
    if False:
        # Wight et al., 1977 method
        x = signal.lfilter([1.,3.41421356,4.87100924,3.41421356,1.], 1.0, x)
        x = x[::2]
        #w, h = signal.freqz([1.,3.41421356,4.87100924,3.41421356,1.])
        #plt.semilogy(w/np.pi, abs(h), 'b')
    
    # Number of segments
    Ns = np.int(np.floor(x.shape[0]/n))
    if x.shape[0] != Ns*n:
        print('Excluding last %d point(s) from calculation.' % (x.shape[0]-Ns*n))
        
    # xr has rows of segment, columns of x
    xr = np.reshape(x[0:Ns*n], (Ns, n))
    sr = np.mean(np.std(xr, axis=1))
    print('%.2f' % sr)
    #print('%d segments of length %d' % (xr.shape[0], xr.shape[1]))

    h = np.hanning(n)
    #h = h/np.mean(h)
    xr = np.multiply(xr, h)

    # Compute FFT of each row
    # In Wight et al. 1977 method, the fft is not used here. Instead, for
    # each segment, the DFT is computed for the 6th and 8th harmonic using 
    # by multiplying each row element by the elements in a pre-computed array
    # of cos/sin terms and then summing. This is done to reduce memory and
    # computation needed in a hardware implementation.
    xr_dft = fft(xr, n) 
    
    # Average across segments. Result has shape of (segment, x(freq))
    xr_dft_ave = np.mean(xr_dft, axis=0) 
    
    f = fftfreq(n)/np.double(factor)
    
    return (f, xr_dft_ave, x)