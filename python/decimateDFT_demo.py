from matplotlib import pyplot as plt
import numpy as np
from numpy.fft import fft, ifft, fftfreq
from decimateDFT import decimateDFT

Nd = 32     # Number of points per decimation
Nt = 32*Nd  # Number of points in time series

tn = 3

if tn == 4:
    x = np.random.randn(Nt)
    xf = fft(x)
    f = fftfreq(Nt)
    I = np.abs(f) > 0.25
    xf[I] = 0
    x = ifft(xf)
    plt.plot(x)
    

if tn == 1:
    t = np.arange(Nt)
    x = np.sin(2.0*np.pi*t/16)
    xpsd = np.abs(fft(x))
    f = fftfreq(Nt)

    plt.clf()
    plt.loglog(f[0:Nt//2],xpsd[0:Nt//2])
    
    (f, x_dft_ave, x_decimated) = decimateDFT(x, Nd, all=True)

    plt.loglog(f[0,0:Nd//2],np.abs(x_dft_ave[0,0:Nd//2]),'.')
    plt.loglog(f[1,0:Nd//2],np.abs(x_dft_ave[1,0:Nd//2]),'.')

if tn == 2:
    
    nx = np.int(Nd/2)
    F = np.array([])
    D = np.array([])
    for d in range(0, f.shape[0]):
        plt.plot(f[d,1:nx], np.abs(x_dft_ave[d,1:nx]))
        F = np.append(F, f[d,6])
        F = np.append(F, f[d,8])
        D = np.append(D, x_dft_ave[d,6])
        D = np.append(D, x_dft_ave[d,8])
        plt.plot(f[d,6], np.abs(x_dft_ave[d,6]),'k*')
        plt.plot(f[d,8], np.abs(x_dft_ave[d,8]),'k*')
        
    I = np.argsort(F)
    F = F[I]
    D = D[I]
    
    if False:
        f = f.flatten()
        x_dft_ave = x_dft_ave.flatten()
        I = np.argsort(f)
        f = f[I]
        x_dft_ave = x_dft_ave[I]
        
if tn == 3:

    Nr = 100

    # Use to get shapes of f    
    (f, x_dft_ave, x_decimated) = decimateDFT(x, Nd, all=True)
    
    X = np.zeros((f.shape[0], f.shape[1], Nr), dtype=complex)
    Y = np.zeros((Nt, Nr), dtype=complex)
    fY = fftfreq(Nt)
    
    for i in range(Nr):
        x = np.random.randn(Nt)
        Y[:,i] = fft(x)
        (f, X[:,:,i], x_decimated) = decimateDFT(x, Nd, all=True)
        
    X = np.abs(X)
    Y = np.abs(Y)
    Xa = np.mean(X, axis=2) # Ave across 
    Ya = np.mean(Y, axis=1)
    plt.clf()
    plt.plot(fY, Ya, '.')
    s = (Nt/Nd)#/np.mean(np.hanning(Nd))
    #for i in range(0, f.shape[1]//2):
    i = 1
    plt.plot(f[:,i], s*Xa[:,i], '*')
    
    plt.legend(np.arange(Xa.shape[1]))


