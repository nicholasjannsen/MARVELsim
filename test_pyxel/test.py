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
pipeline.charge_generation.load_charge.arguments.filename   = 'science.fits'
pipeline.charge_generation.load_charge.arguments.time_scale = 5
pipeline.charge_generation.tars.arguments.seed              = np.random.randint(1e9)

results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)


  # # # photon -> photon
  # # photon_generation:

  #   - name: load_image
  #     func: pyxel.models.photon_generation.load_image
  #     enabled: true
  #     arguments:
  #       image_file: bias.fits      
  #       position: [5000, 5000]
  #       align: "center"
  #       convert_to_photons: true
  #       time_scale: 10

  #   - name: shot_noise
  #     func: pyxel.models.photon_generation.shot_noise
  #     enabled: true
  #     arguments:
  #       type: "poisson"
