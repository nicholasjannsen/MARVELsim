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
  $ python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 -o </path/to/outdir>
  $ python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --rv 5.5 -o </path/to/outdir> 
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
from utilities import errorcode

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
star_group.add_argument('--mag',  metavar='VMAG',   type=str, help='Johnson-Cousin V passband magnitude')
star_group.add_argument('--teff', metavar='KELVIN', type=str, help='Stellar effective temperature [K]')
star_group.add_argument('--logg', metavar='DEX',    type=str, help='Stellar surface gravity [relative log10]')
star_group.add_argument('--z',    metavar='DEX',    type=str, help='Stellar metallicity [Fe/H]')

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

def add_fitsheader(filename, obsmode, exptime):
    # Load and open file
    hdul = fits.open(filename)
    hdr = hdul[0].header
    # Add headers
    hdr.append(('ORIGIN', 'Instituut voor Sterrenkunde, KU Leuven', 'Institution'), end=True)
    hdr.append(('OBSERVAT', 'LaPalma', 'Observatory name'), end=True)
    hdr.append(('TELESCOP', 'MARVEL', 'Telescope name'), end=True)
    hdr.append(('OBSGEO-X', '5327306.5552', 'Cartesian X [meters] GRS80'), end=True)
    hdr.append(('OBSGEO-Y', '-1718448.6952', 'Cartesian Y [meters] GRS80'), end=True)
    hdr.append(('OBSGEO-Z', '3051947.7715', 'Cartesian Z [meters] GRS80'), end=True)
    hdr.append(('OBSERVER', 'Nicholas Jannsen', 'Observer'), end=True)
    hdr.append(('PROG_ID', '0', 'Programme ID'), end=True)
    hdr.append(('INSTRUME', 'MARVEL', 'Instrument'), end=True)
    hdr.append(('FIBMODE', '', 'Fibre mode'), end=True)
    hdr.append(('CREATOR', 'revision_20201117', 'Version of data acquisition system'), end=True)
    hdr.append(('HDRVERS', '20151026', 'Version of FITS header'), end=True)
    hdr.append(('FILENAME', f'{filename}', 'Origina'), end=True)
    hdr.append(('UNSEQ', '995604', 'Unique sequence number'), end=True)
    hdr.append(('OBSMODE', f'{obsmode}', 'Observing mode'), end=True)
    hdr.append(('IMAGETYP', f'{obsmode}', 'Image type'), end=True)
    hdr.append(('EXPTYPE', f'{obsmode}', 'Exposure type'), end=True)
    hdr.append(('COMMENTS', '', 'Free comments by the observer'), end=True)
    hdr.append(('DATE-OBS', '2021-03-09T17:19:37.034931', 'Start of observation'), end=True)
    hdr.append(('DATE-END', '2021-03-09T17:19:37.034167', 'End of observation'), end=True)
    hdr.append(('DATE-AVG', '2021-03-09T17:19:37.034549', 'Midpoint of observation'), end=True)
    hdr.append(('DATE', '2021-03-09T17:20:26.418524', 'Time of file creation'), end=True)
    hdr.append(('EXPTIME', f'{exptime}', 'Exposure time'), end=True)
    hdr.append(('OBJECT', 'CALIBRATION', 'Object name'), end=True)
    hdr.append(('OBJ_RA', '0.0', '[deg] Object RA'), end=True)
    hdr.append(('OBJ_DEC', '0.0', '[deg] Object DEC'), end=True)
    hdr.append(('EQUINOX', '2000.0', 'Equinox of coordinates'), end=True)
    hdr.append(('RADECSYS', 'FK5', 'Coordinate system'), end=True)
    hdr.append(('TEMPM1', '99999.99000000001', '[C] Telescope temperature M1'), end=True)
    hdr.append(('TEMPM1MC', '99999.99000000001', '[C] Telescope temperature mirror cell M1'), end=True)
    hdr.append(('TEMPM2', '99999.99000000001', '[C] Telescope temperature M2'), end=True)
    hdr.append(('TEMPM2E', '7.6', '[C] Telescope temperature M2E'), end=True)
    hdr.append(('TEMPTT', '8.9', '[C] Telescope tube top temperature'), end=True)
    hdr.append(('TEMPTB', '8.5', '[C] Telescope tube center temperature'), end=True)
    hdr.append(('TEMPT032', '99999.99000000001', '[C] Temperature inside REM rack'), end=True)
    hdr.append(('TEMPT106', '99999.99000000001', '[C] Temperature inside RPM rack'), end=True)
    hdr.append(('TEMPT107', '7.9', '[C] Temperature of air in Maia side Nasmyth'), end=True)
    hdr.append(('TEMPT108', '99999.99000000001', '[C] Temperature of Hermes adapter'), end=True)
    hdr.append(('TEMPT033', '99999.99000000001', '[C] Temperature at top op fork, Hermes side'), end=True)
    hdr.append(('TEMPT109', '9.6', '[C] Temperature of air at top of tube'), end=True)
    hdr.append(('TEMPT110', '99999.99000000001', '[C] Temperature of air inside tube'), end=True)
    hdr.append(('TEMPT114', '99999.99000000001', '[C] Dewpoint at top of tube'), end=True)
    hdr.append(('HUMT111', '29.9', '[%] Rel humidity of air at top of tube'), end=True)
    hdr.append(('PCIFILE', '', 'PCI card setup file'), end=True)
    hdr.append(('TIMFILE', '/home/mocs/mocs/config/mocs/marvel/tim-MARVEL-20130205-sp_idle.lod', ''), end=True)
    hdr.append(('UTILFILE', '', 'Utility board setup file'), end=True)
    hdr.append(('DETMODE', 'LGN', 'Controller readout speed/gain setting'), end=True)
    hdr.append(('READMODE', 'L', 'Detector readout mode'), end=True)
    hdr.append(('DETGAIN', '1.2', '[e-/ADU] Detector gain'), end=True)
    hdr.append(('DETBIAS', '2180.0', '[ADU] Expected bias level'), end=True)
    hdr.append(('BINX', '1', 'Binning factor in x'), end=True)
    hdr.append(('BINY', '1', 'Binning factor in y'), end=True)
    hdr.append(('DTM1_1', '1', 'Binning factor in x'), end=True)
    hdr.append(('DTM1_2', '1', 'Binning factor in y'), end=True)
    hdr.append(('WINDOWED', 'FALSE', 'Has the detector been windowed?'), end=True)
    hdr.append(('TEMP_MET', '6.4', '[C]  Temperature Meteo Station'), end=True)
    hdr.append(('HUM_MET', '12.0', '[%]  Rel humidity Meteo Station'), end=True)
    hdr.append(('PRES_MET', '772.0', '[mbar]  Atm pressure Meteo Station'), end=True)
    hdr.append(('WINDAVG', '5.9', '[m/s]  Avg wind speed Meteo Station'), end=True)
    hdr.append(('WINDMAX', '6.9', '[m/s]  Gust wind speed Meteo Station'), end=True)
    hdr.append(('WINDDIR', '27.0', '[deg]  Avg wind direction Meteo Station'), end=True)
    hdr.append(('INSTDATE', '20180724', 'Last intervention in instrument'), end=True)
    hdr.append(('CTRDATE', '20091112', 'Last change in detector controller setup'), end=True)
    hdr.append(('CTR_ID', 'UNKNOWN', 'Controller serial number'), end=True)
    hdr.append(('PCI_ID', 'SN381', 'PCI card serial number'), end=True)
    hdr.append(('TIM_ID', 'UNKNOWN', 'Timing board serial number'), end=True)
    hdr.append(('UTIL_ID', 'UNKNOWN', 'Utility board serial number'), end=True)
    hdr.append(('DETNAME', 'Marvel-Science-GC2', 'Detector name'), end=True)
    hdr.append(('DETTYPE', 'E2V42-90', 'Detector type'), end=True)
    hdr.append(('DETID', 'DET06', 'Detector ID'), end=True)
    hdr.append(('TEMPCCD', '160.0', '[K]  Temperature of detector'), end=True)
    hdr.append(('TEMPCRYO', '83.694', '[K]  Temperature of coldhead'), end=True)
    hdr.append(('PRESH044', '779.232', '[mbar] Pressure in outer room'), end=True)
    hdr.append(('PRESH095', '782.3', '[mbar] Pressure in outer room West'), end=True)
    hdr.append(('HUMH071', '42.355154', '[%] Rel humidity in outer room'), end=True)
    hdr.append(('HUMH072', '35.580974', '[%] Rel humidity on table'), end=True)
    hdr.append(('TEMPH039', '18.008', '[C] Temperature in inner room'), end=True)
    hdr.append(('TEMPH040', '13.838', '[C] Temperature in outer room'), end=True)
    hdr.append(('TEMPH047', '16.996', '[C] HERMES temperature air.camera'), end=True)
    hdr.append(('TEMPH048', '17.851', '[C] HERMES temperature table.center'), end=True)
    hdr.append(('TEMPH050', '17.920', '[C] HERMES temperature grating.mount.top'), end=True)
    hdr.append(('TEMPH051', '17.819', '[C] HERMES temperature fiberexit.mount'), end=True)
    hdr.append(('TEMPH061', '17.859', '[C] HERMES temperature maincoll.glass.top'), end=True)
    hdr.append(('TEMPH052', '18.030', '[C] HERMES temperature maincoll.mount.top'), end=True)
    hdr.append(('TEMPH053', '17.977', '[C] HERMES temperature maincoll.mount.bot'), end=True)
    hdr.append(('TEMPH054', '17.670', '[C] HERMES temperature camera.top.center'), end=True)
    hdr.append(('TEMPH055', '16.483', '[C] HERMES temperature cryostat.front'), end=True)
    hdr.append(('TEMPH056', '16.025', '[C] HERMES temperature cryostat.rear'), end=True)
    hdr.append(('TEMPH057', '17.859', '[C] HERMES temperature grating.glass.center'), end=True)
    hdr.append(('TEMPH058', '17.966', '[C] HERMES temperature maincoll.glass.top'), end=True)
    hdr.append(('TEMPH059', '17.979', '[C] HERMES temperature maincoll.glass.bot'), end=True)
    hdr.append(('BSCALE', '1', ''), end=True)
    hdr.append(('BZERO', '32768', ''), end=True)
    # Write new file with header
    fits.writeto(filename, hdul[0].data, hdr, overwrite=True)

    
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
    ncosmics = int(rate * exptime + np.random.randint(300) * rate)
    pipeline.charge_generation.tars.arguments.particle_number = ncosmics


#==============================================================#
#                      PYECHELLE + PYXEL                       #
#==============================================================#

# Hard-coded parameters
bias_level   = 20398  # [e-] -> 2170 ADU
read_noise   = 52     # [e-] -> 5.5 ADU
exptime_flat = 5
exptime_thar = 5
exptime_thne = 5

# Snippet command for pyechelle
run_marvel = f"pyechelle -s MARVEL_2021_11_22 --fiber 1-5 --bias {bias_level} --read_noise {read_noise}"

# Run with GPUs:
if args.cpu: run_marvel = run_marvel + f" --max_cpu {args.cpu}"

# Run with CUDA requiring NVIDIA cores:
if args.cuda: run_marvel = run_marvel + f" --cuda"

# Create an instance of the pyxel class from the MARVEL specific inputfile
config = pyxel.load("inputfile_marvel.yaml")
exposure = config.exposure
detector = config.ccd_detector
pipeline = config.pipeline

# Set output directory for pyxel
# NOTE the folder "pyxel_dir" cannot be avoided at the moment..
# and nor is it possible to rename the pyxel output files..
output_dir = str(exposure.outputs.output_dir)
pyxel_dir  = output_dir.split('/')[-1]
exposure.outputs.output_dir = pathlib.Path(args.outdir + pyxel_dir)

#------------------------------------#
#            RUN CALIBRATION         #
#------------------------------------#

if args.calibs:  # TODO how many exposures do we need of each calibs?

    # Generate bias (pyechelle)

    for i in range(1,11):
        errorcode('message', '\nSimulating bias')
        # Run pyechelle
        filename_bias = f'{args.outdir}bias_'+f'{i}'.zfill(4)+'.fits'
        command_bias  = f"pyechelle -s MARVEL_2021_11_22 --bias 2000 --read_noise 5.5 --sources Constant -t 0 -o {filename_bias}"
        os.system(command_bias)
        add_fitsheader(filename_bias, 'BIAS', 0)

        # TODO can we do it faster with Pyxel?
        # NOTE We here generate a bias from a shortened dark exposure
        # This is done by multiplying the dark rate with the bias exposure time
        # pipeline.charge_generation.load_charge.enabled = False
        # pipeline.charge_generation.tars.enabled = False
        # pipeline.charge_generation.dark_current.arguments.dark_rate *= float(exptime_bais)
        # pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
    
    # Generate a flat

    for i in range(1,4):
        errorcode('message', '\nSimulating spectral flat')
        # Run pyechelle
        filename_flat = f'{args.outdir}flat_'+f'{i}'.zfill(4)+'.fits'
        command_flat  = run_marvel + f" --sources Constant -t {exptime_flat} -o {filename_flat}"
        os.system(command_flat)
        add_fitsheader(filename_bias, 'FLAT', exptime_flat)
        # Run pyxel
        enable_cosmics(exptime_flat)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_flat
        pipeline.charge_generation.load_charge.arguments.time_scale = float(exptime_flat)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

    # Generate a ThAr arc

    for i in range(1,4):
        errorcode('message', '\nSimulating ThAr arc')
        # Run pyechelle
        filename_thar = f'{args.outdir}thar_'+f'{i}'.zfill(4)+'.fits'
        command_thar  = run_marvel + f" --sources ThAr -t {exptime_thar} -o {filename_thar}"
        os.system(command_thar)
        add_fitsheader(filename_bias, 'THAR', exptime_thar)
        # Run pyxel
        enable_cosmics(exptime_thar)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_thar
        pipeline.charge_generation.load_charge.arguments.time_scale = float(exptime_thar)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

    # Generate a ThNe arc

    for i in range(1,4):
        errorcode('message', '\nSimulating ThNe arc')
        # Run pyechelle
        filename_thne = f'{args.outdir}thne_'+f'{i}'.zfill(4)+'.fits'
        command_thne = run_marvel + f" --sources ThNe -t {exptime_thne} -o {filename_thne}"
        os.system(command_thne)
        add_fitsheader(filename_bias, 'THNE', exptime_thne)
        # Run pyxel
        enable_cosmics(exptime_thne)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_thne
        pipeline.charge_generation.load_charge.arguments.time_scale = float(exptime_thne)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

    # Generate a Etalon & ThAr

    for i in range(1,4):
        errorcode('message', '\nSimulating Etalon & ThAr')
        # Run pyechelle
        filename_wave = f'{args.outdir}wave_'+f'{i}'.zfill(4)+'.fits'
        command_wave  = run_marvel + f" --sources Etalon ThAr ThAr ThAr ThAr --etalon_d=6 -t {exptime_thar} -o {filename_wave}"
        os.system(command_wave)
        add_fitsheader(filename_wave, 'WAVE', exptime_thar)
        # Run pyxel
        enable_cosmics(exptime_thar)
        pipeline.charge_generation.load_charge.arguments.filename   = filename_wave
        pipeline.charge_generation.load_charge.arguments.time_scale = float(exptime_thar)
        pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

#------------------------------------#
#           RUN RV SEQUENCE          #
#------------------------------------#

else:

    # Generate a science frame

    errorcode('message', 'Simulating stellar spectrum\n')
    # Run pyechelle
    if args.index is None: args.index = ''
    if args.rv is None: args.rv = 0
    filename_science = f'{args.outdir}science_'+f'{args.index}'.zfill(4)+'.fits'
    command_science = (f" --sources Phoenix Phoenix Phoenix Phoenix ThAr --etalon_d=6 --d_primary 0.8 --d_secondary 0.1" +
                       f" --phoenix_t_eff {args.teff} --phoenix_log_g {args.logg} --phoenix_z {args.z} --phoenix_alpha 0.0" +
                       f" --phoenix_magnitude {args.mag} --rv {args.rv} --fiber 1-5 -t {args.time} -o {filename_science}")
    os.system(run_marvel + command_science)
    add_fitsheader(filename_science, 'SCIENCE', exptime_thar)
    # Run pyxel
    pipeline.charge_generation.load_charge.arguments.filename   = filename_science
    pipeline.charge_generation.load_charge.arguments.time_scale = float(args.time)
    pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

#==============================================================#
#                          PROLOGUE                            #
#==============================================================#

# Final execution time

toc = datetime.datetime.now()
print(f"\nMARVEL simulations took {toc - tic}")
