#!/usr/bin/env python3

import pyxel
import numpy as np

# Load specific inputfile for MARVEL
config = pyxel.load("inputfile.yaml")

# Setup configurations
exposure = config.exposure        # class Observation
detector = config.ccd_detector    # class CCD
pipeline = config.pipeline        # class DetectionPipeline

# Run pyxel in exposure mode
pipeline.photon_generation.cosmix.arguments.seed                 = np.random.randint(1e9)
pipeline.photon_generation.cosmix.arguments.particles_per_second = 0
pipeline.charge_generation.load_charge.arguments.filename        = 'zeros.fits'
pipeline.charge_generation.load_charge.arguments.time_scale      = 1
pipeline.charge_generation.dark_current.arguments.dark_rate      = 1
#exposure = '[1]'

print(exposure)
exit()



results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
