import numpy as np
from numpy.fft import fft, fftfreq
from astroML.time_series import generate_power_law
from astroML.fourier import PSD_continuous
from matplotlib import pyplot as plt
from decimateDFT import decimateDFT

N = 1024
factor = 100
beta = 1.0
dt = 1.0
t = np.arange(N)
Ns = 1000 # Number of segments

PSD = np.ones((513, Ns))
psd = np.ones((513, Ns))
for i in range(Ns):
    random_state = np.random.RandomState(i)

    # Generate the light curve and compute the PSD
    x = factor * generate_power_law(N, dt, beta, random_state=random_state)

    (fd, x_dft_ave, x_decimated) = decimateDFT(x, 32, all=True)

    f, PSD[:,i] = PSD_continuous(t, x)
    tmp = fft(x)
    psd[:,i] = np.abs(tmp[0:psd.shape[0]])
    psd[:,i] = psd[:,i]*psd[:,i]
    
PSDa = np.mean(PSD, axis=1)
psda = np.mean(psd, axis=1)

#plt.loglog(f[1:], )
#plt.plot(f[1:], PSD[1:,0]/psd[1:,0])
f2 = fftfreq(N)
plt.plot(f[1:], f[1:]*PSDa[1:])
plt.plot(f[1:], 2*f[1:]*psda[1:])
#plt.loglog(f[1:], PSDa[1:])
#plt.loglog(f[1:], psda[1:])
