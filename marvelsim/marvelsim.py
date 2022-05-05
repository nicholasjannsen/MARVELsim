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
  $ python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --rv 5.5 -o </path/to/outdir>
  $ python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --data <rv_data.txt> --cuda -o </path/to/outdir> 
"""

import os
import yaml
import pyxel
import datetime
import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from utilities import errorcode, add_fitsheader
from pathlib import Path
from zipfile import ZipFile

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
    def __init__(self):
        """
        Constructor of the class.
        """

        # HARD-CODED PARAMETERS
        
        self.bias_level   = 2000   # [ADU]
        self.read_noise   = 5      # [ADU] RMS
        
        # PARSED ARGUMENTS
        
        # Input paths
        self.cwd = Path(__file__).parent.resolve()

        # Output paths
        if args.outdir == '.':
            self.outdir = Path.cwd()
        else:
            self.outdir = Path(args.outdir).resolve()

        # Default stellar parameters of a Sun-like star
        if args.mag   is None: args.mag   = 10.
        if args.teff  is None: args.teff  = 5800
        if args.logg  is None: args.logg  = 4.5
        if args.z     is None: args.z     = 0.
        if args.alpha is None: args.alpha = 0.
        if args.rv    is None: args.alpha = 0.

        # Control number of calibration images
        if args.nbias is None: args.nbias = 10
        if args.nflat is None: args.nflat = 5
        if args.nthar is None: args.nthar = 5
        #if args.nthne is None: args.nthne = 5
        if args.nwave is None: args.nwave = 5

        # Control exposure time of calibration images [s]
        if args.tdark is None: args.tdark = 300
        if args.tflat is None: args.tflat = 5
        if args.tthar is None: args.tthar = 30
        #if args.tthne is None: args.tthne = 30
        if args.twave is None: args.twave = 5
        

        
    def init_pyechelle(self):
        """
        Module to initialise PyEchelle
        """
        # Snippet command for pyechelle
        self.run_marvel = f'pyechelle -s MARVEL_2021_11_22 --fiber 1-5 --bias {self.bias_level} --read_noise {self.read_noise}' 

        # Run with normal CPU cores
        if args.cpu:
            self.run_marvel = self.run_marvel + f" --max_cpu {args.cpu}"

        # Run with CUDA requiring NVIDIA hardware and GPUs
        if args.cuda:
            self.run_marvel = self.run_marvel + f" --cuda"

        # Check for RV inputfile
        if args.data:
            try:
                data = np.loadtxt(args.data, delimiter=',', skiprows=1)
            except:
                errorcode('error', r'File do not exist: {args.data}')
            else:
                self.t  = data[:,1]
                args.rv = data[:,2]
        elif args.rv:
            args.rv = [float(args.rv)]
        else:
            args.rv = [0]


    
    def init_pyxel(self):
        """
        Module to initialise Pyxel.
        """        
        
        # Create an instance of the pyxel class from the MARVEL specific inputfile
        filename_inputfile = Path.joinpath(self.cwd, "../inputfiles/inputfile_marvel.yaml")
        config = pyxel.load(filename_inputfile)
        self.exposure = config.exposure
        self.detector = config.ccd_detector
        self.pipeline = config.pipeline

        # Set output directory for pyxel
        # NOTE the folder "pyxel_dir" cannot be avoided at the moment..
        # and nor is it possible to rename the pyxel output files..
        output_dir = str(self.exposure.outputs.output_dir)
        self.pyxel_dir  = output_dir.split('/')[-1]
        self.pyxel_path = Path.joinpath(self.outdir, self.pyxel_dir)
        self.pyxel_file = Path.joinpath(self.pyxel_path, '/detector_image_array_1.fits')
        self.exposure.outputs.output_dir = Path.joinpath(self.outdir, self.pyxel_dir)

        # Finito!
        return self.pyxel_path

    
    
    def enable_cosmics(self, exptime):
        """
        Module to draw cosmic rays from a Poisson distibution scaled to the exposure time.
        """
        # Make sure cosmics are being added
        self.pipeline.photon_generation.cosmix.enabled = True
        # Use spacecraft model for cosmis
        filename_cosmix = Path.joinpath(self.cwd, '../inputfiles/proton_L2_solarMax_11mm_Shielding.txt')
        self.pipeline.photon_generation.cosmix.arguments.spectrum_file = filename_cosmix
        # Set random seed for cosmic rays
        self.pipeline.photon_generation.cosmix.arguments.seed = np.random.randint(1e9)
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
        self.pipeline.photon_generation.cosmix.arguments.particles_per_second = ncosmics

        

    def fetch_nimg(self, imgtype):
        """
        Module fetch the number of exposure of image type.
        """
        if imgtype == 'bias':    nimg = args.nbias
        if imgtype == 'flat':    nimg = args.nflat
        if imgtype == 'thar':    nimg = args.nthar
        if imgtype == 'wave':    nimg = args.nwave
        if imgtype == 'science': nimg = len(args.rv)
        # Finito!
        return nimg


    
    def fetch_exptime(self, imgtype):
        """
        Module fetch the exposure time of image type.
        """
        if imgtype == 'bias':    exptime = 0
        if imgtype == 'flat':    exptime = args.tflat
        if imgtype == 'thar':    exptime = args.tthar
        if imgtype == 'wave':    exptime = args.twave
        if imgtype == 'science': exptime = args.time
        # Finito!
        return exptime

    

    def compress_data(self, filename, filepath):
        # Compress file
        print(f'Compressing {filename}')
        filepath = str(filepath) 
        with ZipFile(f'{filepath[:-5]}.zip', 'w') as zipfile:
            zipfile.write(filepath, f'{filename}')
        # Remove uncompressed file
        os.remove(filepath)
            
    #--------------------------------------------#
    #                  PYECHELLE                 #  
    #--------------------------------------------#

    def cmd_pyechelle(self, imgtype, filepath, i):

        if imgtype == 'bias':
            cmd = (f'pyechelle -s MARVEL_2021_11_22 --sources Constant -t 0' +
                   f' --bias {self.bias_level} --read_noise {self.read_noise} -o {filepath}')
        if imgtype == 'flat':
            cmd = (self.run_marvel +
                   f' --sources Constant --constant_intensity 0.01' +
                   f' -t {args.tflat} -o {filepath}')
        if imgtype == 'thar':
            cmd = (self.run_marvel +
                   f' --sources ThAr -t {args.tthar} -o {filepath}')
        # if imgtype == 'thne':
        #     cmd = (run_marvel +
        #            f' --sources ThNe -t {args.tthne} -o {filepath}')
        if imgtype == 'wave':
            cmd = (self.run_marvel +
                   f' --sources Etalon ThAr ThAr ThAr ThAr --etalon_d=6' + 
                   f' -t {args.tthar} -o {filepath}')
        if imgtype == 'science':
            cmd = (self.run_marvel +
                   ' --etalon_d=6 --d_primary 0.8 --d_secondary 0.1' +
                   f' --sources Phoenix Phoenix Phoenix Phoenix ThAr' +
                   f' --phoenix_t_eff {args.teff} --phoenix_log_g {args.logg}' +
                   f' --phoenix_z {args.z} --phoenix_alpha {args.alpha}' +
                   f' --phoenix_magnitude {args.mag}' +
                   f' --rv {args.rv[i]} -t {args.time} -o {filepath}')
        # Finito!
        return cmd


    
    def run_pyechelle(self, imgtype, fitstype):
        """
        Module to generate spectra with PyEchelle.
        """        
        for i in range(1, self.fetch_nimg(imgtype)+1):
            errorcode('message', f'\nSimulating {imgtype} with PyEchelle')
            # Run pyechelle
            filename = f'{imgtype}_'+f'{i}'.zfill(4)+'.fits'
            filepath = Path.joinpath(self.outdir, filename)
            os.system(self.cmd_pyechelle(imgtype, filepath, i-1))

            # TODO can we do it faster with Pyxel?
            # NOTE We here generate a bias from a shortened dark exposure
            # This is done by multiplying the dark rate with the bias exposure time
            # pipeline.charge_generation.load_charge.enabled = False
            # pipeline.charge_generation.tars.enabled = False
            # pipeline.charge_generation.dark_current.arguments.dark_rate *= float(exptime_bais)
            # pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)
            if imgtype == 'bias':
                add_fitsheader(filepath, fitstype, 0)
                # Compress file
                if args.zip:
                    self.compress_data(filename, filepath)

    #--------------------------------------------#
    #                    PYXEL                   #  
    #--------------------------------------------#

    def run_pyxel(self, imgtype, fitstype):
        """
        Module to add CCD effects spectra with Pyxel.
        """
        for i in range(1, self.fetch_nimg(imgtype)+1):
            errorcode('message', f'\nSimulating {imgtype} with Pyxel')
            # Fetch filenames.
            filename = f'{imgtype}_'+f'{i}'.zfill(4)+'.fits'
            filepath = Path.joinpath(self.outdir, filename)
            # Run pyxel
            self.enable_cosmics(self.fetch_exptime(imgtype))
            self.pipeline.charge_generation.load_charge.arguments.filename   = filepath
            self.pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(args.time)
            pyxel.exposure_mode(exposure=self.exposure, detector=self.detector, pipeline=self.pipeline)
            # Swap files
            os.remove(filepath)
            os.system(f'mv {self.pyxel_file} {filepath}')
            # Lastly add header
            print('Adding fits-header')
            add_fitsheader(filepath, fitstype, args.time)
            # Compress file
            if args.zip:
                self.compress_data(filename, filepath)


                
    def run_pyxel_cpu(self, imgtype, fitstype):
        """
        Module to add CCD effects spectra with Pyxel.
        """
        errorcode('message', f'\nSimulating {imgtype} with Pyxel')
        # Fetch filenames.
        filename = f'{imgtype}_'+f'{args.dex}'.zfill(4)+'.fits'
        filepath = Path.joinpath(self.outdir, filename)
        # Run pyxel
        self.enable_cosmics(self.fetch_exptime(imgtype))
        self.pipeline.charge_generation.load_charge.arguments.filename   = filepath
        self.pipeline.charge_generation.load_charge.arguments.time_scale = 5 #float(args.time)
        pyxel.exposure_mode(exposure=self.exposure, detector=self.detector, pipeline=self.pipeline)
        # Swap files
        os.remove(filepath)
        os.system(f'mv {self.pyxel_file} {filepath}')
        # Lastly add header
        print('Adding fits-header')
        add_fitsheader(filepath, fitstype, args.time)
        # Compress file
        if args.zip:
            self.compress_data(filename, filepath)

                
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
cal_group.add_argument('--nthar', metavar='NUM', type=int, help='Number of ThAr exposures (default: 5)')
cal_group.add_argument('--nflat', metavar='NUM', type=int, help='Number of Flat exposures (default: 5)')
cal_group.add_argument('--nwave', metavar='NUM', type=int, help='Number of Etalon + ThAr exposures  (default: 5)')
cal_group.add_argument('--tdark', metavar='NUM', type=int, help='Exposure time of Dark (default: 300 s)')
cal_group.add_argument('--tflat', metavar='NUM', type=int, help='Exposure time of Flat (default: 5 s)')
cal_group.add_argument('--tthar', metavar='NUM', type=int, help='Exposure time of ThAr (default: 30 s)')
cal_group.add_argument('--twave', metavar='NUM', type=int, help='Exposure time of Etalon + ThAr (default: 30 s)')

hpc_group = parser.add_argument_group('PERFORMANCE')
hpc_group.add_argument('--data', metavar='PATH', type=str, help='Path to include RV file')
hpc_group.add_argument('--dex',  metavar='NAME', type=int, help='Index for running Pyxel on CPUs')
hpc_group.add_argument('--cpu',  metavar='INT',  type=str, help='Maximum number of CPU cores used order-wise parallel computing')
hpc_group.add_argument('--cuda', action='store_true', help='Flag to use CUDA NVIDIA hardware for raytracing (makes cpu flag obsolete')
hpc_group.add_argument('--zip',  action='store_true', help='Flag to zip output files.')

args = parser.parse_args()

# Create instance of class
m = marvelsim()
m.init_pyechelle()

if args.calibs:
    # Run pyechelle
    m.run_pyechelle('bias', 'BIAS')
    m.run_pyechelle('flat', 'FLAT')
    m.run_pyechelle('thar', 'THAR')
    m.run_pyechelle('wave', 'WAVE')
    # Run Pyxel
    pyxel_path = m.init_pyxel()
    m.run_pyxel('flat', 'FLAT')
    m.run_pyxel('thar', 'THAR')
    m.run_pyxel('wave', 'WAVE')
    os.rmdir(pyxel_path)
    
else:

    if args.cuda or args.cpu:
        # Run pyechelle alone with either CUDA or CPUs 
        m.run_pyechelle('science', 'SCIENCE')
    
    if args.dex:
        # Run pyxel alone with CPUs
        pyxel_path = m.init_pyxel()
        m.run_pyxel_cpu('science', 'SCIENCE')
        os.rmdir(pyxel_path)

    else:
        # Run pyechelle
        m.run_pyechelle('science', 'SCIENCE')
        # Run pyxel
        pyxel_path = m.init_pyxel()
        m.run_pyxel('science', 'SCIENCE')
        os.rmdir(pyxel_path)

# Final execution time
toc = datetime.datetime.now()
print(f"\nMARVEL simulations took {toc - tic}")
