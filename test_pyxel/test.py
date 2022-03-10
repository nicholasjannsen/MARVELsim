#!/usr/bin/env python3

import pyxel

# Load specific inputfile for MARVEL
config = pyxel.load("inputfile.yaml")

# Setup configurations
exposure = config.exposure        # class Observation
detector = config.ccd_detector    # class CCD
pipeline = config.pipeline        # class DetectionPipeline

# Run pyxel in exposure mode
#print(exposure.outputs.output_dir)
#print(pipeline.charge_generation.load_charge.arguments)
print(pipeline.photon_generation.cosmix)
exit()
pipeline.charge_generation.load_charge.arguments.filename   = 'science.fits'
#pipeline.charge_generation.load_charge.arguments.time_scale = 300
#pipeline.charge_generation.tars.arguments.seed              = np.random.randint(1e9)

results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

print(results)
