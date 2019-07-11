#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 22:35:16 2019

@author: robertweigel
"""

import numpy as np
from numpy.fft import ifft, fft, fftfreq
from scipy import signal
from matplotlib import pyplot as plt

y = np.array([1, 1+1j,1+1j,1+1j,1+1j,1+1j])
yft = ifft(y)