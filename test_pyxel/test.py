#!/usr/bin/env python3

import pyxel
import numpy as np
import matplotlib.pyplot as plt

# Load specific inputfile for MARVEL
config = pyxel.load("test.yaml")

# Setup configurations
exposure = config.exposure        # class Observation
detector = config.ccd_detector    # class CCD
pipeline = config.pipeline        # class DetectionPipeline

# Run pyxel in exposure mode
pipeline.charge_generation.load_charge.arguments.filename   = 'flat_0001.fits'
pipeline.charge_generation.load_charge.arguments.time_scale = 5
pipeline.charge_generation.tars.arguments.seed              = np.random.randint(1e9)
results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
