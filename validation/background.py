#!/usr/bin/env python3

"""
This simulation tool is designed specifically for the upcoming MARVEL
spectroscopic survey consisting of 4 small robotic telescopes each of
a 0.8 primary mirror diameter, and combined by fibers to a highly stable
high resolution spectrographs. This simulator uses the PyEchelle code
to simulate realistic MARVEL spectra using a ray tracing technique, and
the Pyxel code that adds various important CCD effects. Seen above this
software can be used to simulate standard calibrated spectra and stellar
spectra including a RV signal. Used in combination with High Performance
Computing (HPC) this makes it easy to simulate a time series of spectra.
"""

# Defaults
import os

# External
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import biweight_location, mad_std
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from scipy import interpolate

#==============================================================#
#                           UTILITIES                          #
#==============================================================#

# Load data
hdul = fits.open('00906451_HRF_FF.fits')
hdr  = hdul[0].header
data = hdul[0].data

# Reduce some statistics
#print(f'Median   value: {np.median(data)}')
#print(f'Biweight value: {biweight_location(data)}')
#print(f'MAD std  value: {mad_std(data)}')

# Perform 
sigma_clip = SigmaClip(sigma=3.0)
bkg_estimator = MedianBackground()
bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

# Statistics for background image
print(f'Bg median: {bkg.background_median}')
print(f'Bg rms   : {bkg.background_rms_median}')

# New image
bg    = bkg.background
data1 = data - bg

# Plot data
norm = ImageNormalize(stretch=SqrtStretch())
fig, ax = plt.subplots(1, 3, figsize=(15,10))
ax[0].imshow(data,  norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
ax[1].imshow(bg,    norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
ax[2].imshow(data1, norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
plt.show()


# Interpolate image
x = np.linspace(0, 1, bg.shape[0])
y = np.linspace(0, 1, bg.shape[1])
f = interpolate.interp2d(y, x, bg, kind='cubic')
x1 = np.linspace(0, 1, 10560)
y1 = np.linspace(0, 1, 10560)
bg1 = f(y1, x1)

# Plot data
plt.figure()
plt.imshow(bg1, norm=norm, origin='lower', cmap='Greys_r', interpolation='nearest')
plt.show()

