>#!/usr/bin/env python3

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
from utilities import errorcode, add_fitsheader

# Turn off warnings
import warnings
warnings.filterwarnings("ignore")

# Monitor script speed
tic = datetime.datetime.now()

#==============================================================#
#                      PYECHELLE + PYXEL                       #
#==============================================================#

class marvelsim(object):
    """
    This class has the purpose
    """
    
    # INITILIZE THE CLASS:
    def __init__(self, args):
        """
        Constructor of the class.
        """

        # HARD-CODED PARAMETERS
        
        self.bias_level   = 2000   # [ADU]
        self.read_noise   = 5      # [ADU] RMS
        
        # PARSED ARGUMENTS
        
        # Make it possible to use bash syntax for PWD
        if args.outdir == '.':
            self.outdir = os.getcwd()
        else:
            self.outdir = args.outdir

        # Control number of calibration images
        if args.nbias is None: self.nbias = 10
        else: self.nbias = args.nbias
        if args.ndark is None: self.ndark = 10
        else: self.ndark = args.ndark
        if args.nflat is None: self.flat = 5
        else: self.nflat = args.nflat
        if args.nthar is None: self.nthar = 5
        else: self.nthar = args.nthar
        #if args.nthne is None: self.nthne = 5
        #else: self.nthne = args.nthne
        if args.nwave is None: self.nwave = 5
        else: self.nwave = args.nwave

        # Control exposure time of calibration images [s]
        if args.tdark is None: self.tdark = 300
        else: self.ndark = args.ndark
        if args.tflat is None: self.tflat = 5
        else: self.tflat = args.tflat
        if args.tthar is None: self.tthar = 30
        else: self.tthar = args.tthar
        #if args.tthne is None: self.tthne = 30
        #else: self.thne = args.tthne
        if args.twave is None: self.twave = 5
        else: self.twave = args.twave

        
        


        

    def init_pyechelle(self, args):
        """
        Module to initialise PyEchelle
        """
        # Snippet command for pyechelle
        self.run_marvel = f'pyechelle -s MARVEL_2021_11_22 --fiber 1-5 --bias {self.bias_level} --read_noise {self.read_noise}' 

        # Run with GPUs:
        if args.cpu:
            self.run_marvel = self.run_marvel + f" --max_cpu {args.cpu}"

        # Run with CUDA requiring NVIDIA cores:
        if args.cuda:
            self.run_marvel = self.run_marvel + f" --cuda"



    def init_pyxel(self, args):
        """
        Module to initialise Pyxel.
        """
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



    def enable_cosmics(self, exptime):
        """
        Module to draw a random number distribution of cosmics scaled to the exposure time.
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


        

    def run_calibs_pyechelle(self, args):
        """
        Module to run generate calibration data with PyEchelle.
        """
        
        # Generate bias (pyechelle)

        for i in range(1,self.nbias):
            errorcode('message', '\nSimulating bias')
            # Run pyechelle
            filename_bias = f'{args.outdir}/bias_'+f'{i}'.zfill(4)+'.fits'
            command_bias  = (f'pyechelle -s MARVEL_2021_11_22 --sources Constant -t 0' +
                             f' --bias {self.bias_level} --read_noise {self.read_noise} -o {filename_bias}')
            os.system(command_bias)
            add_fitsheader(filename_bias, 'BIAS', 0)

        # Generate a ThAr arc

        for i in range(1,self.nthar):
            errorcode('message', '\nSimulating ThAr arc')
            # Run pyechelle
            filename_thar = f'{args.outdir}/thar_'+f'{i}'.zfill(4)+'.fits'
            command_thar  = self.run_marvel + f" --sources ThAr -t {self.exptime_thar} -o {filename_thar}"
            os.system(command_thar)

        # Generate a ThNe arc

        # for i in num_calibs:
        #     errorcode('message', '\nSimulating ThNe arc')
        #     # Run pyechelle
        #     filename_thne = f'{args.outdir}/thne_'+f'{i}'.zfill(4)+'.fits'
        #     command_thne = run_marvel + f" --sources ThNe -t {self.exptime_thne} -o {filename_thne}"
        #     os.system(command_thne)

        # Generate a Etalon & ThAr

        for i in range(1,self.nwave):
            errorcode('message', '\nSimulating Etalon & ThAr')
            # Run pyechelle
            filename_wave = f'{args.outdir}/wave_'+f'{i}'.zfill(4)+'.fits'
            command_wave  = (self.run_marvel + f' --sources Etalon ThAr ThAr ThAr ThAr --etalon_d=6 ' + 
                             '-t {self.exptime_thar} -o {filename_wave}')
            os.system(command_wave)

        # Generate a flat

        for i in range(1,self.nflat):
            errorcode('message', '\nSimulating spectral flat')
            # Run pyechelle
            filename_flat = f'{args.outdir}/flat_'+f'{i}'.zfill(4)+'.fits'
            command_flat  = self.run_marvel + f" --sources Constant --constant_intensity 0.01 -t {self.exptime_flat} -o {filename_flat}"
            os.system(command_flat)
        


            

    def run_calibs_pyxel(self):
        """
        Module to add CCDs effects to the calibrated data produced by PyEchelle.
        """

        # TODO can we do it faster with Pyxel?
        # NOTE We here generate a bias from a shortened dark exposure
        # This is done by multiplying the dark rate with the bias exposure time
        # pipeline.charge_generation.load_charge.enabled = False
        # pipeline.charge_generation.tars.enabled = False
        # pipeline.charge_generation.dark_current.arguments.dark_rate *= float(exptime_bais)
        # pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
        
        # Generate a ThAr arc
        errorcode('message', '\n[pyxel]: Simulating ThAr arc')
        for i in range(1,self.nbias):
            # Run pyxel
            self.enable_cosmics(self.exptime_thar)
            filename_thar = f'{args.outdir}/thar_'+f'{i}'.zfill(4)+'.fits'
            pipeline.charge_generation.load_charge.arguments.filename   = filename_thar
            pipeline.charge_generation.load_charge.arguments.time_scale = 5.0 #float(exptime_thar)
            pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
            # Swap files
            os.remove(filename_thar)
            os.system(f'mv {pyxel_file} {filename_thar}')
            # Add header
            add_fitsheader(filename_thar, 'THAR', exptime_thar)

        # Generate a ThNe arc

        # for i in num_calibs:
        #     errorcode('message', '\nSimulating ThNe arc')
        #     # Run pyxel
        #     enable_cosmics(exptime_thne)
        #     filename_thne = f'{args.outdir}/thne_'+f'{i}'.zfill(4)+'.fits'
        #     pipeline.charge_generation.load_charge.arguments.filename   = filename_thne
        #     pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(exptime_thne)
        #     pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
        #     # Swap files
        #     os.remove(filename_thne)
        #     os.system(f'mv {pyxel_file} {filename_thne}')
        #     # Add header
        #     add_fitsheader(filename_thne, 'THNE', exptime_thne)

        # Generate a Etalon & ThAr

        for i in range(1,self.nthar):
            errorcode('message', '\nSimulating Etalon & ThAr')
            # Run pyxel
            self.enable_cosmics(exptime_thar)
            filename_wave = f'{args.outdir}/wave_'+f'{i}'.zfill(4)+'.fits'
            pipeline.charge_generation.load_charge.arguments.filename   = filename_wave
            pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(exptime_wave)
            pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
            # Swap files
            os.remove(filename_wave)
            os.system(f'mv {pyxel_file} {filename_wave}')
            # Add header
            add_fitsheader(filename_wave, 'WAVE', exptime_wave)

        # Generate a flat

        for i in range(1,self.nflat):
            errorcode('message', '\nSimulating spectral flat')
            # Run pyxel
            self.enable_cosmics(exptime_flat)
            filename_flat = f'{args.outdir}/flat_'+f'{i}'.zfill(4)+'.fits'
            pipeline.charge_generation.load_charge.arguments.filename   = filename_flat
            pipeline.charge_generation.load_charge.arguments.time_scale = 5.0 #float(exptime_flat)
            pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
            # Swap files
            os.remove(filename_flat)
            os.system(f'mv {pyxel_file} {filename_flat}')
            # Add header
            add_fitsheader(filename_flat, 'FLAT', exptime_flat)




            
    def run_science_pyechelle(self, args):
        """
        Module to run PyEchelle for star spectra.
        """
        
        # Check if all stellar parameters are present
        star = [args.teff, args.logg, args.z, args.alpha]
        if None in star:
            errorcode('error', 'One or more star parameters are missing!')
        
        # Generate a science frame
        
        errorcode('message', '\nSimulating stellar spectrum with PyEchelle\n')
        # Run pyechelle
        if args.rv is None: args.rv = 0
        if args.index is None: args.index = 1
        filename_science = f'{args.outdir}/science_'+f'{args.index}'.zfill(4)+'.fits'
        command_science = (f" --sources Phoenix Phoenix Phoenix Phoenix ThAr --etalon_d=6 --d_primary 0.8 --d_secondary 0.1" +
                           f" --phoenix_t_eff {args.teff} --phoenix_log_g {args.logg} --phoenix_z {args.z} --phoenix_alpha {args.alpha}" +
                           f" --phoenix_magnitude {args.mag} --rv {args.rv} -t {args.time} -o {filename_science}")
        os.system(self.run_marvel + command_science)



        

    def run_science_pyxel(self, args):
        """
        Module to run PyEchelle for star spectra.
        """
                
        errorcode('message', '\nSimulating stellar spectrum with PyEchelle\n')
        # Run pyxel
        self.enable_cosmics(args.time)
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
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

software = '\nThe MARVEL Spectroscopy Simulator -> MARVELsim'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

parser.add_argument('-o', '--outdir', metavar='PATH', type=str, help='Output directory to store simulations')

obs_group = parser.add_argument_group('OBSERVATION')
obs_group.add_argument('-t', '--time', metavar='SEC', type=float, help='Exposure time of stellar observation [s]')
obs_group.add_argument('--mag',   metavar='VMAG',     type=str,   help='Johnson-Cousin V passband magnitude')
obs_group.add_argument('--teff',  metavar='KELVIN',   type=str,   help='Stellar effective temperature [K]')
obs_group.add_argument('--logg',  metavar='DEX',      type=str,   help='Stellar surface gravity [relative log10]')
obs_group.add_argument('--z',     metavar='DEX',      type=str,   help='Stellar metallicity [Fe/H]')
obs_group.add_argument('--alpha', metavar='DEX',      type=str,   help='Stellar Alpha element abundance [alpha/H]')
obs_group.add_argument('--rv',    metavar='M/S',      type=str,   help='Radial Velocity shift of star due to exoplanet [m/s]')

cal_group = parser.add_argument_group('CALIBRATION')
cal_group.add_argument('-c', '--calibs', action='store_true', help='Flag to simulate a calibration dataset (default: False)')
cal_group.add_argument('--nbias', metavar='NUM', type=int, help='Number of Bias exposures (default: 10)')
cal_group.add_argument('--ndark', metavar='NUM', type=int, help='Number of Dark exposures (default: 10)')
cal_group.add_argument('--nthar', metavar='NUM', type=int, help='Number of ThAr exposures (default:  5)')
cal_group.add_argument('--nflat', metavar='NUM', type=int, help='Number of Flat exposures (default:  5)')
cal_group.add_argument('--nwave', metavar='NUM', type=int, help='Numexp of Etalon + ThAr  (default:  5)')

cal_group.add_argument('--tdark', metavar='NUM', type=int, help='Exposure time of Dark (default: 300 s)')
cal_group.add_argument('--tflat', metavar='NUM', type=int, help='Exposure time of Flat (default:   5 s)')
cal_group.add_argument('--tthar', metavar='NUM', type=int, help='Exposure time of ThAr (default:  30 s)')
cal_group.add_argument('--twave', metavar='NUM', type=int, help='Exptime Etalon+ThAr   (default:  30 s)')

hpc_group = parser.add_argument_group('PERFORMANCE')
hpc_group.add_argument('--index', metavar='INT', type=str, help='Integer index used for parallel computations')
hpc_group.add_argument('--cpu',   metavar='INT', type=str, help='Maximum number of CPU cores used order-wise parallel computing')
hpc_group.add_argument('--cuda',  action='store_true', help='NVIDIA hardware using CUDA used for raytracing (makes cpu flag obsolete')
hpc_group.add_argument('--hpc',   action='store_true', help='Flag to tell software when running on the hpc')

args = parser.parse_args()

# Run software
m = marvelsim(args)

# Calibration
if args.calibs:
    if args.cuda and args.hpc:
        m.init_pyechelle(args)
        m.run_calibs_pyechelle(args)
    elif args.hpc:
        m.init_pyxel(args)
        m.run_calibs_pyxel(args)
    else:
        m.init_pyechelle(args)
        m.run_calibs_pyechelle(args)
        m.init_pyxel(args)
        m.run_calibs_pyxel(args)

# Science spectra
else:
    if args.cuda and args.hpc:
        m.init_pyechelle(args)
        m.run_science_pyechelle(args)
    elif args.hpc:
        m.init_pyxel(args)
        m.run_science_pyxel(args)
    else:
        m.init_pyechelle(args)
        m.run_science_pyechelle(args)
        m.init_pyxel(args)
        m.run_science_pyxel(args)

# Remove pyxel output folder
if args.cuda is None:
    os.rmdir(args.outdir + '/' + pyxel_dir)

# Final execution time
toc = datetime.datetime.now()
print(f"\nMARVEL simulations took {toc - tic}")
