#!/usr/bin/env python3

import pyxel
import numpy as np

# Load specific inputfile for MARVEL
config = pyxel.load("../../inputfiles/inputfile_marvel.yaml")
#config = pyxel.load("inputfile.yaml")

# Setup configurations
exposure = config.exposure        # class Observation
detector = config.ccd_detector    # class CCD
pipeline = config.pipeline        # class DetectionPipeline

# Exposure class
pipeline.photon_generation.cosmix.arguments.spectrum_file          = "../../inputfiles/proton_L2_solarMax_11mm_Shielding.txt"
pipeline.photon_generation.cosmix.arguments.seed                   = np.random.randint(1e9)
pipeline.photon_generation.cosmix.arguments.particles_per_second   = 0
pipeline.charge_generation.load_charge.arguments.filename          = 'bias_slow.fits'
pipeline.charge_generation.simple_dark_current.arguments.dark_rate = 0.01
pipeline.charge_generation.load_charge.arguments.time_scale        = 1
exposure.readout.times  = [1]

# Run pyxel in exposure mode
results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
