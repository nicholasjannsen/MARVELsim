#!/usr/bin/env python3

"""
Small script to generate the background.
"""

# Defaults
import os

# External
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.background import Background2D, MedianBackground
from scipy.interpolate import RectBivariateSpline

#==============================================================#
#                           UTILITIES                          #
#==============================================================#

# Load data
hdul = fits.open('00906451_HRF_FF.fits')
hdr  = hdul[0].header
data = hdul[0].data

# Reduce some statistics
mean, median, std = sigma_clipped_stats(data, sigma=3.0)
print('Statistics BEFORE background subtracted')
print(f'mean   value: {mean}')
print(f'median value: {median}')
print(f'std    value: {std}')

# Determine the 2D background image
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

# Statistics for background image
print('Statistics AFTER background subtracted')
print(f'Bg median: {bkg.background_median}')
print(f'Bg rms   : {bkg.background_rms_median}')

# New image
bg    = bkg.background
data1 = data - bg
norm = ImageNormalize(stretch=SqrtStretch())

# Plot data
fig, ax = plt.subplots(1, 3, figsize=(15,10))
ax[0].imshow(data,  norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
ax[1].imshow(bg,    norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
ax[2].imshow(data1, norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
plt.show()

# Interpolate image to MARVEL size
x = np.linspace(0, 1, bg.shape[0])
y = np.linspace(0, 1, bg.shape[1])
f = RectBivariateSpline(x, y, bg) #, kind='cubic')
x1 = np.linspace(0, 1, 10560)
y1 = np.linspace(0, 1, 10560)
bg1 = f(x1, y1)

# Plot data
plt.figure()
plt.imshow(bg1, norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
plt.show()

