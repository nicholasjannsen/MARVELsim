#!/usr/bin/env python3

import pyxel
import numpy as np
import matplotlib.pyplot as plt

def scale_linear(inputArray, scale_min=None, scale_max=None):
    """
    Performs linear scaling of the input np array.

    PARAMETERS
    ----------
	inputArray : np array
	    Image data array
	scale_min : float
	    Minimum data value
	scale_max : float
	    Maximum data value

    RETURN
    ------
	imageData : np array
        Image data array
	"""
    # Avoid data is being overwritten
    imageData = np.array(inputArray, copy=True)

    # Select min and max image values as default
    if scale_min is None: scale_min = imageData.min()
    if scale_max is None: scale_max = imageData.max()

    # Scale data array and return
    imageData = imageData.clip(min=scale_min, max=scale_max)
    imageData = (imageData - scale_min) / (scale_max - scale_min)
    indices = np.where(imageData < 0)
    imageData[indices] = 0.0
    indices = np.where(imageData > 1)
    imageData[indices] = 1.0
    return imageData


def FITS(img, sigma=2, cmap='Blues_r', colorbar=True):
    """
    Automatical scale image and plot 2D array.
    Nothing returned but a plt.show() is needed to show the image.

    PARAMETERS
    ----------
	img : np arrap
	    Image data array
	sigma : int, float
	    Scaling parameter as too the number of std sigma's is used.
	cmap : str
       Colormap from matplotlib.pyplot library.
    colorbar : bool
       Whether or not a colorbar should be plotted along.
	"""

    # Find min and max scale for image
    img_min = img.mean() - sigma*img.std()
    img_max = img.mean() + sigma*img.std()

    plt.imshow(scale_linear(img, img_min, img_max), cmap=cmap, origin='lower')

    if colorbar:
        # Prepare for scaled colorbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        ax = plt.gca()
        im = ax.imshow(scale_linear(img, img_min, img_max), cmap=cmap, origin='lower')

        # create an axes on the right side of ax. The width of cax will be 5%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax  = divider.append_axes('right', size='5%', pad=0.05)  # right, left, top, bottom

        # Make colorbar and append ylabel and axis labels:
        cbar = plt.colorbar(im, cax=cax)# orientation='horizontal')
        cbar.ax.set_ylabel('Normalized Counts')

#-------------------

# Load specific inputfile for MARVEL
config = pyxel.load("test.yaml")

# Setup configurations
exposure = config.exposure        # class Observation
detector = config.ccd_detector    # class CCD
pipeline = config.pipeline        # class DetectionPipeline

# Run pyxel in exposure mode
results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

# Show image
# image = results.image.values[0]
# plt.figure(figsize=(12,12))
# FITS(image, cmap='hot')
# plt.show()
