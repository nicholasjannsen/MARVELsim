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
  $ python simulator-marvel.py --calibs -o </path/to/outdir>
  $ python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --rv 5.5 -o </path/to/outdir>
  $ python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --rv $rv --index $i -o </path/to/outdir> 
"""

import os
import pyxel
import pathlib
import datetime
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from utilities import errorcode, add_fitsheader

# Turn off warnings
import warnings
warnings.filterwarnings("ignore")

# Monitor script speed
tic = datetime.datetime.now()

#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

software = '\nThe MARVEL Spectroscopy Simulator -> MARVELsim'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

parser.add_argument('-c', '--calibs', action='store_true', help='Flag to invoke a set of calibration data. Default: False')
parser.add_argument('-o', '--outdir', metavar='OUTDIR', type=str, help='Output directory to store simulations')

obs_group = parser.add_argument_group('OBSERVATION')
obs_group.add_argument('-t', '--time', metavar='SEC', type=float, help='Exposure time of stellar observation [s]')

star_group = parser.add_argument_group('STAR')
star_group.add_argument('--mag',   metavar='VMAG',   type=str, help='Johnson-Cousin V passband magnitude')
star_group.add_argument('--teff',  metavar='KELVIN', type=str, help='Stellar effective temperature [K]')
star_group.add_argument('--logg',  metavar='DEX',    type=str, help='Stellar surface gravity [relative log10]')
star_group.add_argument('--z',     metavar='DEX',    type=str, help='Stellar metallicity [Fe/H]')
star_group.add_argument('--alpha', metavar='DEX',    type=str, help='Stellar Alpha element abundance [alpha/H]')

planet_group = parser.add_argument_group('EXOPLANET')
planet_group.add_argument('--rv', metavar='M/S', type=str, help='Radial Velocity shift of star due to exoplanet [m/s]')

hpc_group = parser.add_argument_group('PERFORMANCE')
hpc_group.add_argument('--index', metavar='INT', type=str, help='Integer index used for parallel computations')
hpc_group.add_argument('--cpu',   metavar='INT', type=str, help='Maximum number of CPU cores used order-wise parallel computing')
hpc_group.add_argument('--cuda',  action='store_true', help='NVIDIA hardware using CUDA used for raytracing (makes cpu flag obsolete')

args = parser.parse_args()

#==============================================================#
#                           UTILITIES                          #
#==============================================================#
    
def enable_cosmics(exptime):
    """
    Draw a random number distribution of cosmics scaled to the exposure time.
    """
    # Make sure cosmics are being added
    pipeline.charge_generation.tars.enabled = True
    # Set random seed for cosmic rays
    pipeline.charge_generation.tars.arguments.seed = np.random.randint(1e9)
    # Benchmark 100 cosmics to an exposure time of 300 seconds
    #--------- testing
    #r = 100/300  # Rate
    #k = np.random.randint(100)
    #l = r * exptime
    #Ncosmics = l**k * np.exp(-l) / np.math.factorial(k)
    #print(Ncosmics)
    #exit()
    rate = 100/300.
    ncosmics = int(rate * exptime + np.random.randint(500) * rate)
    pipeline.charge_generation.tars.arguments.particle_number = ncosmics


#==============================================================#
#                      PYECHELLE + PYXEL                       #
#==============================================================#

# Hard-coded parameters
bias_level   = 2000   # [ADU]
read_noise   = 5      # [ADU] RMS
exptime_thar = 30     # [s]
exptime_thne = 30     # [s]
exptime_wave = 10     # [s]
exptime_flat = 5      # [s]
num_calibs = range(1,6)

# Make it possible to use bash syntax for PWD
if args.outdir == '.': args.outdir = os.getcwd()

# Snippet command for pyechelle
run_marvel = f'pyechelle -s MARVEL_2021_11_22 --fiber 1-5 --bias {bias_level} --read_noise {read_noise}' 

# Run with GPUs:
if args.cpu: run_marvel = run_marvel + f" --max_cpu {args.cpu}"

# Run with CUDA requiring NVIDIA cores:
if args.cuda: run_marvel = run_marvel + f" --cuda"

# Create an instance of the pyxel class from the MARVEL specific inputfile
config = pyxel.load("inputfiles/inputfile_marvel.yaml")
exposure = config.exposure
detector = config.ccd_detector
pipeline = config.pipeline

# Set output directory for pyxel
# NOTE the folder "pyxel_dir" cannot be avoided at the moment..
# and nor is it possible to rename the pyxel output files..
output_dir = str(exposure.outputs.output_dir)
pyxel_dir  = output_dir.split('/')[-1]
pyxel_file = args.outdir + '/' + pyxel_dir + '/detector_image_array_1.fits'
exposure.outputs.output_dir = pathlib.Path(args.outdir + '/' + pyxel_dir)

#------------------------------------#
#            RUN CALIBRATION         #
#------------------------------------#

if args.calibs:  # TODO how many exposures do we need of each calibs?

    # Generate bias (pyechelle)

    for i in range(1,11):
        errorcode('message', '\nSimulating bias')
        # Run pyechelle
        filename_bias = f'{args.outdir}bias_'+f'{i}'.zfill(4)+'.fits'
        command_bias  = (f'pyechelle -s MARVEL_2021_11_22 --sources Constant -t 0' +
                         f' --bias {bias_level} --read_noise {read_noise} -o {filename_bias}')
        os.system(command_bias)
        add_fitsheader(filename_bias, 'BIAS', 0)

        # TODO can we do it faster with Pyxel?
        # NOTE We here generate a bias from a shortened dark exposure
        # This is done by multiplying the dark rate with the bias exposure time
        # pipeline.charge_generation.load_charge.enabled = False
        # pipeline.charge_generation.tars.enabled = False
        # pipeline.charge_generation.dark_current.arguments.dark_rate *= float(exptime_bais)
        # pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
    
    # Generate a ThAr arc

    for i in num_calibs:
        errorcode('message', '\nSimulating ThAr arc')
        # Run pyechelle
        filename_thar = f'{args.outdir}thar_'+f'{i}'.zfill(4)+'.fits'
        command_thar  = run_marvel + f" --sources ThAr -t {exptime_thar} -o {filename_thar}"
        os.system(command_thar)
        # Run pyxel
        enable_cosmics(exptime_thar)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_thar
        pipeline.charge_generation.load_charge.arguments.time_scale = 5.0 #float(exptime_thar)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
        # Swap files
        os.remove(filename_thar)
        os.system(f'mv {pyxel_file} {filename_thar}')
        # Add header
        add_fitsheader(filename_thar, 'THAR', exptime_thar)

    # Generate a ThNe arc

    for i in num_calibs:
        errorcode('message', '\nSimulating ThNe arc')
        # Run pyechelle
        filename_thne = f'{args.outdir}thne_'+f'{i}'.zfill(4)+'.fits'
        command_thne = run_marvel + f" --sources ThNe -t {exptime_thne} -o {filename_thne}"
        os.system(command_thne)
        # Run pyxel
        enable_cosmics(exptime_thne)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_thne
        pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(exptime_thne)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
        # Swap files
        os.remove(filename_thne)
        os.system(f'mv {pyxel_file} {filename_thne}')
        # Add header
        add_fitsheader(filename_thne, 'THNE', exptime_thne)

    # Generate a Etalon & ThAr

    for i in num_calibs:
        errorcode('message', '\nSimulating Etalon & ThAr')
        # Run pyechelle
        filename_wave = f'{args.outdir}wave_'+f'{i}'.zfill(4)+'.fits'
        command_wave  = run_marvel + f" --sources Etalon ThAr ThAr ThAr ThAr --etalon_d=6 -t {exptime_thar} -o {filename_wave}"
        os.system(command_wave)
        # Run pyxel
        enable_cosmics(exptime_thar)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_wave
        pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(exptime_wave)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
        # Swap files
        os.remove(filename_wave)
        os.system(f'mv {pyxel_file} {filename_wave}')
        # Add header
        add_fitsheader(filename_wave, 'WAVE', exptime_wave)

    # Generate a flat

    for i in num_calibs:
        errorcode('message', '\nSimulating spectral flat')
        # Run pyechelle
        filename_flat = f'{args.outdir}/flat_'+f'{i}'.zfill(4)+'.fits'
        command_flat  = run_marvel + f" --sources Constant --constant_intensity 0.01 -t {exptime_flat} -o {filename_flat}"
        os.system(command_flat)
        # Run pyxel
        enable_cosmics(exptime_flat)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_flat
        pipeline.charge_generation.load_charge.arguments.time_scale = 5.0 #float(exptime_flat)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
        # Swap files
        os.remove(filename_flat)
        os.system(f'mv {pyxel_file} {filename_flat}')
        # Add header
        add_fitsheader(filename_flat, 'FLAT', exptime_flat)
        
#------------------------------------#
#           RUN RV SEQUENCE          #
#------------------------------------#

else:

    # Generate a science frame

    errorcode('message', '\nSimulating stellar spectrum\n')
    # Run pyechelle
    if args.rv is None: args.rv = 0
    if args.index is None: args.index = 1
    filename_science = f'{args.outdir}/science_'+f'{args.index}'.zfill(4)+'.fits'
    command_science = (f" --sources Phoenix Phoenix Phoenix Phoenix ThAr --etalon_d=6 --d_primary 0.8 --d_secondary 0.1" +
                       f" --phoenix_t_eff {args.teff} --phoenix_log_g {args.logg} --phoenix_z {args.z} --phoenix_alpha {args.alpha}" +
                       f" --phoenix_magnitude {args.mag} --rv {args.rv} -t {args.time} -o {filename_science}")
    os.system(run_marvel + command_science)
    # Run pyxel
    enable_cosmics(args.time)
    pipeline.charge_generation.load_charge.arguments.filename   = filename_science
    pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(args.time)
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
    # Swap files
    os.remove(filename_science)
    os.system(f'mv {pyxel_file} {filename_science}')
    os.rmdir(args.outdir + '/' + pyxel_dir)
    # Lastly add header
    add_fitsheader(filename_science, 'SCIENCE', args.time)

#==============================================================#
#                          PROLOGUE                            #
#==============================================================#

# Remove pyxel output folder
os.rmdir(args.outdir + '/' + pyxel_dir)

# Final execution time
toc = datetime.datetime.now()
print(f"\nMARVEL simulations took {toc - tic}")
