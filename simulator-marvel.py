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

User examples:
  $ python  
  $ python 
  $ python 
"""

import os
import time
import pyxel
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from colorama import Fore, Style

# Turn off warnings
import warnings
warnings.filterwarnings("ignore")

# Monitor script speed
tic = time.time()

#==============================================================#
#                           UTILITIES                          #
#==============================================================#

def errorcode(API, message):
    """
    This function allows to colour code error messages within a code.
    """
    if API == 'software':
        print(Style.BRIGHT + Fore.GREEN + message + Style.RESET_ALL)
    if API == 'message':
        print(Style.BRIGHT + message + Style.RESET_ALL)
    if API == 'warning':
        print(Style.BRIGHT + Fore.YELLOW + '[Warning]: ' + message + Style.RESET_ALL)
    if API == 'error':
        print(Style.BRIGHT + Fore.RED + '[Error]: ' + message + Style.RESET_ALL)
        exit()

#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

software = '\nThe MARVEL Spectroscopy Simulator -> MARVELsim'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

parser.add_argument('outdir', type=str, help='Output directory to store simulations')
parser.add_argument('-c', '--calibs', action='store_true', help='Flag to invoke a set of calibration data. Default: False')

obs_group = parser.add_argument_group('OBSERVATION')
obs_group.add_argument('-t', '--time', action='store_true', help='Exposure time of stellar observation [s]')

star_group = parser.add_argument_group('STAR')
star_group.add_argument('--mag',  action='store_true', help='Johnson-Cousin V passband magnitude')
star_group.add_argument('--teff', action='store_true', help='Stellar effective temperature [K]')
star_group.add_argument('--logg', action='store_true', help='Stellar surface gravity [relative log10]')
star_group.add_argument('--z',    action='store_true', help='Stellar metallicity [Fe/H]')

planet_group = parser.add_argument_group('EXOPLANET')
planet_group.add_argument('--rv', action='store_true', help='Radial Velocity shift of star due to exoplanet [m/s]')

args = parser.parse_args()

#==============================================================#
#                      PYECHELLE + PYXEL                       #
#==============================================================#

# Snippet command for pyechelle

run_marvel = "pyechelle -s MARVEL_2021_11_22"

# Create an instance of the pyxel class from the MARVEL specific inputfile

config = pyxel.load("inputfile_marvel.yaml")
exposure = config.exposure
detector = config.ccd_detector
pipeline = config.pipeline

# Set output directory for pyxel

exposure.outputs.output_dir = args.outdir + 'flat_pyxel.fits'

#------------------------------------#
#            RUN CALIBRATION         #
#------------------------------------#

if args.calibs:  # TODO how many exposures do we need of each calibs?

    # Hard code exposure times
    exptime_flat = 1
    exptime_thar = 5
    exptime_thne = 5
    exptime_bais = 0.001  # NOTE Typical clocking of CCD

    # Generate a flat

    errorcode('message', '\nSimulating spectral flat')
    # Run pyechelle
    command_flat = run_marvel + f" --fiber 2-5 --sources Constant -t {exptime_flat} -o {args.outdir}flat.fits"
    os.system(command_flat)
    # Run pyxel
    pipeline.charge_generation.load_charge.arguments.filename   = f'{args.outdir}flat.fits'
    pipeline.charge_generation.load_charge.arguments.time_scale = exptime_flat  # TODO testing
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
    exit()

    # Generate a ThAr arc

    errorcode('message', '\nSimulating ThAr arc')
    # Run pyechelle
    command_thar = run_marvel + f" --fiber 2-5 --sources ThAr -t {exptime_thar} -o {args.outdir}thar.fits"
    #os.system(command_thar)
    # Run pyxel
    pipeline.charge_generation.load_charge.arguments.filename   = f'{args.outdir}thar.fits'
    pipeline.charge_generation.load_charge.arguments.time_scale = float(exptime_thar)
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

    # Generate a ThNe arc

    errorcode('message', '\nSimulating ThNe')
    # Run pyechelle
    command_thne = run_marvel + f" --fiber 2-5 --sources ThNe -t {exptime_thne} -o {args.outdir}thne.fits"
    os.system(command_thne)
    # Run pyxel
    pipeline.charge_generation.load_charge.arguments.filename   = f'{args.outdir}thne.fits'
    pipeline.charge_generation.load_charge.arguments.time_scale = float(exptime_thne)
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

    # Generate a arc frame
    # TODO needed in real life but they are already included in science frames?

    #errorcode('message', 'Simulating spectral arc\n')
    #command_wave = run_marvel + " --sources Etalon --etalon_d=6 -t 10 -o marvel_etalon.fits"
    #os.syste,(command_wave)

    # Generate bias with pyxel only:
    # NOTE We here generate a bias from a shortened dark exposure
    # This is done by multiplying the dark rate with the bias exposure time

    errorcode('message', '\nSimulating bias')
    pipeline.charge_generation.load_charge.enabled = False
    pipeline.charge_generation.tars.enabled = False
    pipeline.charge_generation.dark_current.arguments.dark_rate *= float(exptime_bais)
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

#------------------------------------#
#           RUN RV SEQUENCE          #
#------------------------------------#

else:

    # Generate a science frame

    errorcode('message', 'Simulating stellar spectrum\n')
    # Run pyechelle
    command_science = (f" --sources Phoenix Phoenix Phoenix Phoenix ThAr --etalon_d=6 --d_primary 0.8 --d_secondary 0.1" +
                       f" --phoenix_t_eff {args.Teff} --phoenix_log_g {args.logg} --phoenix_z {args.Z} --phoenix_alpha 0.0" +
                       f" --phoenix_magnitude {args.Vmag} --fiber 1-5 -t {args.time} -o marvel_science.fits")
    os.system(run_marvel + command_science)
    # Run pyxel
    pipeline.charge_generation.load_charge.arguments.filename   = f'{args.outdir}/thne.fits'
    pipeline.charge_generation.load_charge.arguments.time_scale = float(f'{args.time}')
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

#==============================================================#
#                          PROLOGUE                            #
#==============================================================#

# Final execution time

toc = time.time()
print(f"MARVEL simulations took {toc - tic} s")
