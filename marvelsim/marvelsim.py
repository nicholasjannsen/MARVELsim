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
"""

# Defaults
import os
import glob
import shutil
import zipfile
import warnings
import datetime
import argparse
import subprocess

# External
import yaml
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from pathlib import Path

# PyEchelle
from pyechelle.simulator import Simulator
from pyechelle.sources import Constant, ThAr, ThNe, Etalon, Phoenix
from pyechelle.spectrograph import ZEMAX
from pyechelle.telescope import Telescope
from pyechelle.CCD import CCD

# Pyxel
import pyxel

# MARVELsim
from utilities import errorcode, add_fitsheader

# Turn off warnings
warnings.filterwarnings("ignore")

# Monitor script speed
tic = datetime.datetime.now()

#==============================================================#
#                        MARVELsim CLASS                       #          
#==============================================================#

class marvelsim(object):
    """
    Class to control all functions related to the MARVELsim software.
    """
    
    def __init__(self):
        """
        Constructor of the class.
        """

        # CCD PARAMETERS
        
        # Image dimentions
        self.bit = 65536
        self.dim = 10560

        # Resolution mode
        # TODO only high (H) for now -> future use ultr-high (U)
        self.resmode = 'H'
        
        # Default readout mode
        if args.readmode is None: self.readmode = 'slow'
        else: self.readmode = args.readmode
        
        # Properties for readout modes
        if self.readmode in {'slow', 'slowmax'}:
            self.speed   = 1000         # [kHz]
            self.gain_ac = 3.0          # [e/ADU]
            self.gain_dc = 21.9e-6      # [V/ADU]
        elif self.readmode in {'fast', 'fastmax'}:
            self.speed   = 50           # [kHz]
            self.gain_ac = 9.4          # [e/ADU]
            self.gain_dc = 65.8e-6      # [V/ADU]
        else:
            errorcode('error', f'Invalid readout mode: {args.readmode}')

        # Read noise RMS [e-]
        if   self.readmode == 'fast':    self.read_noise = 2.5  # Nominal fast-readout [e-]
        elif self.readmode == 'fastmax': self.read_noise = 4.0  # Maximum fast-readout [e-]
        elif self.readmode == 'slow':    self.read_noise = 5.0  # Nominal slow-readout [e-]
        elif self.readmode == 'slowmax': self.read_noise = 7.0  # Maximum slow-readout [e-]

        # ADC input range [V]
        self.adc_range = [0, self.gain_dc*self.bit]  

        # Gain of output amplifier [V -> V] 
        self.gain_amp = 5  

        # Charge readout sensitivity [V/e-]
        self.readsen = 7.0e-6
        
        # Bias offset [ADU] and bias level [e-]
        self.bias       = 1000
        self.bias_level = int(self.bias * self.gain_ac)
        
        # Initialise random number generator
        # NOTE if not parsed it use after the clock
        if not args.seed: self.rng = np.random.default_rng()
        else:             self.rng = np.random.default_rng(args.seed)
        self.seed = self.rng.integers(1e9, size=1)[0]

        # PARSED ARGUMENTS
        
        # Input paths
        self.cwd = Path(__file__).parent.resolve()

        # Output paths
        if args.outdir == '.':
            self.outdir = Path.cwd()
        elif args.outdir:
            self.outdir = Path(args.outdir).resolve()
        else:
            errorcode('error', f'Provide a output path! Use -o </path/to/output>')

        # Create output path if it doesn't exist
        self.outdir.mkdir(parents=True, exist_ok=True)
            
        # Make sure not to overwrite rv
        self.rv = args.rv
        
        # Default stellar parameters of a Sun-like star
        if args.time  is None: args.time  = 900.
        if args.mag   is None: args.mag   = [10.]
        if args.teff  is None: args.teff  = [5800]
        if args.logg  is None: args.logg  = [4.5]
        if args.z     is None: args.z     = [0.]
        if args.alpha is None: args.alpha = [0.]
        if args.rv    is None: args.rv    = [0.]

        # Control number of calibration images
        if args.calibs:
            if args.nbias is None: args.nbias = 10
            if args.ndark is None: args.ndark = 5
            if args.nflat is None: args.nflat = 5
            if args.nthar is None: args.nthar = 5
            if args.nwave is None: args.nwave = 5
            #if args.nthne is None: args.nthne = 5
        else:
            if args.nbias is None: args.nbias = 1
            if args.ndark is None: args.ndark = 1
            if args.nflat is None: args.nflat = 1
            if args.nthar is None: args.nthar = 1
            if args.nwave is None: args.nwave = 1

        # Control exposure time of calibration images [s]
        if args.tdark is None: args.tdark = 900
        if args.tflat is None: args.tflat = 1.5
        if args.tthar is None: args.tthar = 15
        if args.twave is None: args.twave = 15
        #if args.tthne is None: args.tthne = 30
        

        
    def fetch_nimg(self, imgtype):
        """
        Module fetch the number of exposure of image type.
        """
        try: self.nexp
        except: self.nexp = 1
        else: pass
        
        if imgtype == 'BBBBB': nimg = args.nbias
        if imgtype == 'DDDDD': nimg = args.ndark
        if imgtype == 'FFFFF': nimg = args.nflat
        if imgtype == 'TTTTT': nimg = args.nthar
        if imgtype == 'ETTTT': nimg = args.nwave
        if imgtype == 'ESSSS': nimg = self.nexp
        if imgtype == 'TSSSS': nimg = self.nexp

        # Finito!
        return nimg


    
    def fetch_exptime(self, imgtype):
        """
        Module fetch the exposure time of image type.
        """
        if imgtype == 'BBBBB': exptime = 0.001
        if imgtype == 'DDDDD': exptime = args.tdark
        if imgtype == 'FFFFF': exptime = args.tflat
        if imgtype == 'TTTTT': exptime = args.tthar
        if imgtype == 'ETTTT': exptime = args.twave
        if imgtype == 'ESSSS': exptime = args.time
        if imgtype == 'TSSSS': exptime = args.time
        # Finito!
        return exptime



    def fetch_filename(self, imgtype, resmode, exptime):
        """
        Module fetch the exposure time of image type.
        """
        now  = datetime.datetime.now()
        date = f'{now.strftime("%Y")}{now.strftime("%m")}{now.strftime("%d")}'
        time = datetime.datetime.now().strftime("%H%M%S")
        filename = f'{date}T{time}_{imgtype}_{resmode}_'+f'{int(exptime)}'.zfill(4)+'.fits'
        # Finito!
        return filename
    
    
    
    def compress_data(self, filename, filepath):
        """
        Module to give full file access and compress a fits file.
        """
        print(f'Compressing {filename}')
        filepath = str(filepath)
        # Open and write zip file
        with zipfile.ZipFile(f'{filepath[:-5]}.zip', 'w', compression=zipfile.ZIP_DEFLATED) as fzip:
            fzip.write(filepath, filename)
            # Give full read/write permission
            os.system(f'chmod 755 {filepath[:-5]}.zip')
        # Remove uncompressed file
        os.remove(filepath)


        
    #--------------------------------------------#
    #                 BIAS IMAGES                #  
    #--------------------------------------------#
        
    def run_bias(self, imgtype):
        """
        Module to generate both bias and dark data arrays.
        Dark arrays are here identical to bias arrays but darks will be
        post-processed later on by Pyxel.
        """
        for i in range(1, self.fetch_nimg(imgtype)+1):
            errorcode('message', f'\nSimulating {imgtype} with Numpy')
            # Setup paths
            exptime  = self.fetch_exptime(imgtype)
            filename = self.fetch_filename(imgtype, self.resmode, exptime)
            filepath = self.outdir / filename
            # Create data array
            bias = self.rng.normal(self.bias_level, self.read_noise, size=(self.dim, self.dim)).astype(int)
            # Save to fits file
            hdul = fits.HDUList([fits.PrimaryHDU(bias)])
            hdul.writeto(filepath)


    #--------------------------------------------#
    #                  PYECHELLE                 #  
    #--------------------------------------------#
            
    def init_pyechelle(self):
        """
        Module to initialise PyEchelle
        """

        # Configure setup for MARVEL
        self.sim = Simulator(ZEMAX('MARVEL_2021_11_22'))
        self.sim.set_ccd(1)
        self.sim.set_fibers([1, 2, 3, 4, 5])
        self.sim.set_telescope(Telescope(0.8, 0.1))
        
        # Add bias level and read noise
        self.sim.set_bias(self.bias_level)
        self.sim.set_read_noise(self.read_noise)
        
        # Run with CUDA requiring NVIDIA hardware and GPUs
        if args.cuda:
            self.sim.set_cuda(True, self.seed)

        # Check for multiple target stars
        if args.science or args.debug=="ESSSS":
            if len(args.mag) == 1: args.mag = np.ones(4) * args.mag
            elif len(args.mag) == 4: pass
            else: errorcode('error', 'mag takes only 1 or 4 values!')

            if len(args.teff) == 1: args.teff = np.ones(4) * args.teff
            elif len(args.teff) == 4: pass
            else: errorcode('error', 'teff takes only 1 or 4 values!')
            # Make sure that Teff is expressed as integers (PyEchelle requirement)
            args.teff = [int(i) for i in args.teff]
            
            if len(args.logg) == 1: args.logg = np.ones(4) * args.logg
            elif len(args.logg) == 4: pass
            else: errorcode('error', 'logg takes only 1 or 4 values!')
            
            if len(args.z) == 1: args.z = np.ones(4) * args.z
            elif len(args.z) == 4: pass
            else: errorcode('error', 'z takes only 1 or 4 values!')

            if len(args.alpha) == 1: args.alpha = np.ones(4) * args.alpha
            elif len(args.alpha) == 4: pass
            else: errorcode('error', 'alpha takes only 1 or 4 values!')

            # Check for RV inputfile
            if args.data:
                # For parallisation
                try: data = np.loadtxt(args.data, delimiter=',', skiprows=1)
                except: errorcode('error', r'File do not exist: {args.data}')
                else: args.rv = data[:,2:]
                # Number of exposures
                self.nexp = len(args.rv)
            else:
                # Multiple targets, single CPU
                if len(args.rv) == 1: args.rv = np.ones(4) * args.rv
                elif len(args.rv) == 4: args.rv = np.array(args.rv)
                else: errorcode('error', 'alpha takes only 1 or 4 values!')
                # Only one exposure
                self.nexp = 1

                
    
    def run_pyechelle(self, imgtype):
        """
        Module to generate spectra with PyEchelle.
        """
        
        # Set exposure times
        exptime = self.fetch_exptime(imgtype)
        self.sim.set_exposure_time(exptime)

        # Adjust etalon intensity
        n_etalon = 3e5
        if not int(exptime) == 0: 
            n_etalon /= exptime
            
        # Run loop over each exposure
        for i in range(self.fetch_nimg(imgtype)):
            errorcode('message', f'\nSimulating {imgtype} with PyEchelle')
            
            # PyEchelle in Python needs to be re-initialized
            #self.init_pyechelle()

            if imgtype == 'FFFFF':
                self.sim.set_sources(Constant(intensity=0.01))
            
            if imgtype == 'TTTTT':
                self.sim.set_sources(ThAr())
                
            if imgtype == 'ETTTT':
                self.sim.set_sources([Etalon(d=6, n_photons=n_etalon), ThAr(), ThAr(), ThAr(), ThAr()])

            if imgtype == 'ESSSS':
                # Set etalon + star(s)
                self.sim.set_sources([Etalon(d=6, n_photons=n_etalon),
                                      Phoenix(t_eff=args.teff[0], log_g=args.logg[0], z=args.z[0], alpha=args.alpha[0], magnitude=args.mag[0]),
                                      Phoenix(t_eff=args.teff[1], log_g=args.logg[1], z=args.z[1], alpha=args.alpha[1], magnitude=args.mag[1]),
                                      Phoenix(t_eff=args.teff[2], log_g=args.logg[2], z=args.z[2], alpha=args.alpha[2], magnitude=args.mag[2]),
                                      Phoenix(t_eff=args.teff[3], log_g=args.logg[3], z=args.z[3], alpha=args.alpha[3], magnitude=args.mag[3])])

            if imgtype == 'TSSSS':
                # Set etalon + star(s)
                self.sim.set_sources([ThAr(),
                                      Phoenix(t_eff=args.teff[0], log_g=args.logg[0], z=args.z[0], alpha=args.alpha[0], magnitude=args.mag[0]),
                                      Phoenix(t_eff=args.teff[1], log_g=args.logg[1], z=args.z[1], alpha=args.alpha[1], magnitude=args.mag[1]),
                                      Phoenix(t_eff=args.teff[2], log_g=args.logg[2], z=args.z[2], alpha=args.alpha[2], magnitude=args.mag[2]),
                                      Phoenix(t_eff=args.teff[3], log_g=args.logg[3], z=args.z[3], alpha=args.alpha[3], magnitude=args.mag[3])])

                
                # Activate atmospheric transmission for target(s)
                self.sim.set_atmospheres([False, True, True, True, True])

                # Set radial velocity of stellar target(s) + Etalon
                if np.ndim(args.rv) == 1:
                    # Select value for one target
                    args.rv = [0., args.rv[0], args.rv[1], args.rv[2], args.rv[3]]
                else:
                    # Dataset: Check for 1 common RV signal or 4 different RVs
                    if args.rv.shape[1]==1:   args.rv = [0., args.rv[i-1],    args.rv[i-1],    args.rv[i-1],    args.rv[i-1]]
                    elif args.rv.shape[1]==4: args.rv = [0., args.rv[0][i-1], args.rv[1][i-1], args.rv[2][i-1], args.rv[3][i-1]]
                    else: errorcode('error', 'RV data file need 1 or 4 columns of RVs!')
                self.sim.set_radial_velocities(args.rv)

            # Set outputfile location
            filename = self.fetch_filename(imgtype, self.resmode, exptime)
            filepath = self.outdir / filename
            self.sim.set_output(filepath, overwrite=True)

            # Run simulation
            self.sim.run()

            
    #--------------------------------------------#
    #                    PYXEL                   #  
    #--------------------------------------------#

    def init_pyxel(self):
        """
        Module to initialise Pyxel.
        """
        # We copy and alter the YAML file to correct the "output_folder" location
        # Define orginal YAML file and the copy output YAML file
        ifile = self.cwd / "../inputfiles/inputfile_marvel.yaml"
        ofile = self.outdir / "inputfile_marvel.yaml"
        # First copy YAML to avoid overwriting the original file
        shutil.copy2(ifile, ofile)
        # Load the data within the YAML file
        stream = open(str(ofile), 'r')
        data   = yaml.full_load(stream)
        # Alter the output file location
        data['exposure']['outputs']['output_folder'] = str(self.outdir)
        # Overwrite YAML file
        with open(str(ofile), 'w') as yaml_file:
            yaml_file.write(yaml.dump(data, default_flow_style=False))
            
        # Create an instance of the pyxel class from the MARVEL specific inputfile
        config = pyxel.load(str(ofile))
        self.exposure = config.exposure
        self.detector = config.ccd_detector
        self.pipeline = config.pipeline

        # Remove YAML input file again (pathlib syntax)
        ofile.unlink()
        
        # Set output directory for pyxel
        # NOTE the folder "pyxel_dir" cannot be avoided at the moment..
        # and nor is it possible to rename the pyxel output files..
        # This is our workaround:
        output_dir = str(self.exposure.outputs.output_dir)
        self.pyxel_dir  = output_dir.split('/')[-1]
        self.pyxel_path = self.outdir / self.pyxel_dir
        self.pyxel_file = self.pyxel_path / 'detector_image_array_1.fits'
        self.exposure.outputs.output_dir = self.pyxel_path

        # Set CCD characteristics from class
        self.detector.characteristics.adc_voltage_range = self.adc_range
        self.detector.characteristics.pre_amplification = self.gain_amp
        
        # Finito!
        return self.pyxel_path

    
    
    def enable_cosmics(self, imgtype):
        """
        Module to draw cosmic rays from a Poisson distibution scaled to the exposure time.
        """
        # Make sure cosmics are being added
        if imgtype == 'BBBBB':
            self.pipeline.photon_generation.cosmix.enabled = False
        else:
            self.pipeline.photon_generation.cosmix.enabled = True

        # Use spacecraft model for cosmis
        filename_cosmix = self.cwd / '../inputfiles/proton_L2_solarMax_11mm_Shielding.txt'
        self.pipeline.photon_generation.cosmix.arguments.spectrum_file = filename_cosmix

        # Set random seed for cosmic rays
        self.pipeline.photon_generation.cosmix.arguments.seed = self.seed

        # Benchmark 100 cosmics to an exposure time of 300 seconds
        self.pipeline.photon_generation.cosmix.arguments.particles_per_second = 1.
        #--------- testing
        #r = 100/300  # Rate
        #k = np.random.randint(100)
        #l = r * exptime
        #Ncosmics = l**k * np.exp(-l) / np.math.factorial(k)
        #ncosmics = int(rate * exptime + np.random.randint(500) * rate)



    
    def run_pyxel(self, imgtype):
        """
        Module to add CCD effects spectra with Pyxel.
        """

        # Create a new folder XXXXX to seperate files
        imgdir = self.outdir / imgtype
        nimg   = self.fetch_nimg(imgtype)
        if nimg > 0: imgdir.mkdir(parents=True, exist_ok=True)

        # Loop over each image
        for i in range(nimg):
            errorcode('message', f'\nSimulating {imgtype} with Pyxel')

            # Set exposure time and time scale (to be identical)
            # NOTE Pyxel do not handle 0s -> Set correct time in case user sets 0s
            exptime = self.fetch_exptime(imgtype)
            if int(exptime) == 0: exptime = self.fetch_exptime('BBBBB')
            self.exposure.readout.times = [exptime]
            self.pipeline.charge_generation.load_charge.arguments.time_scale = exptime            

            # Fetch and set filename
            filepath = Path(glob.glob(str(self.outdir) + f'/*{imgtype}*.fits')[i])
            self.pipeline.charge_generation.load_charge.arguments.filename = filepath

            # Enable cosmic rays
            self.enable_cosmics(imgtype)
            
            # Run pyxel
            pyxel.exposure_mode(exposure=self.exposure, detector=self.detector, pipeline=self.pipeline)

            # Remove pyEchelle file not debug
            if args.debug is None:
                os.remove(filepath)

            # Rename Pyxel-PyEchelle file to the old PyEchelle file name
            filename = filepath.name
            filepath = imgdir / filename
            os.system(f'mv {self.pyxel_file} {filepath}')

            # Add fits header
            print('Adding fits-header')
            add_fitsheader(filepath, filename, imgtype, exptime,
                           self.readmode, self.speed, self.gain_ac, self.bias, 
                           self.readsen, self.gain_amp, self.adc_range)
            # Give full read/write permission
            os.system(f'chmod 755 {filepath}')

            # Compress file
            if args.zip:
                self.compress_data(filename, filepath)


                
    def run_pyxel_cpu(self, imgtype):
        """
        Module to add CCD effects spectra with Pyxel.
        """
        errorcode('message', f'\nSimulating {imgtype} with Pyxel')
        # Fetch filenames.
        filename = f'{imgtype}_'+f'{args.dex}'.zfill(4)+'.fits'
        filepath = self.outdir / filename
        # Run pyxel
        exptime = self.fetch_exptime(imgtype)
        self.enable_cosmics(imgtype)
        self.pipeline.charge_generation.load_charge.arguments.filename   = filepath
        self.pipeline.charge_generation.load_charge.arguments.time_scale = exptime
        self.exposure.readout.times = [exptime]
        pyxel.exposure_mode(exposure=self.exposure, detector=self.detector, pipeline=self.pipeline)
        # Swap files
        os.remove(filepath)
        os.system(f'mv {self.pyxel_file} {filepath}')
        # Add fits header
        print('Adding fits-header')
        add_fitsheader(filepath, imgtype, exptime, self.readmode, self.bias, self.gain_ac, self.speed)
        # Give full read/write permission
        os.system(f'chmod 755 {filepath}')
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
obs_group.add_argument('-t', '--time',  metavar='SEC',   type=float, help='Exposure time of stellar observation [s]')
obs_group.add_argument('-f', '--frame', metavar='XXXXX', type=str,   help='Frame to be simulated (e.g. ETTTT)')
obs_group.add_argument('-d', '--debug', metavar='XXXXX', type=str,   help='Debug mode saving also the PyEchelle frame (e.g. ETTTT)')

sci_group = parser.add_argument_group('SCIENCE MODE')
sci_group.add_argument('-s', '--science', action='store_true', help='Flag to simulate stellar spectra')
sci_group.add_argument('--mag',   nargs='*', type=float, help='Johnson-Cousin V passband magnitude (max 4 values)')
sci_group.add_argument('--teff',  nargs='*', type=int,   help='Stellar effective temperature [K] (max 4 values)')
sci_group.add_argument('--logg',  nargs='*', type=float, help='Stellar surface gravity [relative log10] (max 4 values)')
sci_group.add_argument('--z',     nargs='*', type=float, help='Stellar metallicity [Fe/H] (max 4 values)')
sci_group.add_argument('--alpha', nargs='*', type=float, help='Stellar Alpha element abundance [alpha/H] (max 4 values)')
sci_group.add_argument('--rv',    nargs='*', type=float, help='Radial Velocity shift of star(s) due to exoplanet [m/s] (max 4 values)')

cal_group = parser.add_argument_group('CALIBRATION MODE')
cal_group.add_argument('-c', '--calibs', action='store_true', help='Flag to simulate calibration data')
cal_group.add_argument('--nbias', metavar='NUM', type=int, help='Number of Bias exposures (default: 10)')
cal_group.add_argument('--ndark', metavar='NUM', type=int, help='Number of Dark exposures (default: 10)')
cal_group.add_argument('--nthar', metavar='NUM', type=int, help='Number of ThAr exposures (default:  5)')
cal_group.add_argument('--nflat', metavar='NUM', type=int, help='Number of Flat exposures (default:  5)')
cal_group.add_argument('--nwave', metavar='NUM', type=int, help='Number of Etalon + ThAr  (default:  5)')
cal_group.add_argument('--tdark', metavar='SEC', type=float, help='Exposure time of Dark (default: 300 s)')
cal_group.add_argument('--tflat', metavar='SEC', type=float, help='Exposure time of Flat (default:   5 s)')
cal_group.add_argument('--tthar', metavar='SEC', type=float, help='Exposure time of ThAr (default:  30 s)')
cal_group.add_argument('--twave', metavar='SEC', type=float, help='Exptime Etalon + ThAr (default:  30 s)')

ccd_group = parser.add_argument_group('CCD DETECTOR')
ccd_group.add_argument('--readmode', metavar='NAME', type=str, help='Readout mode: fast, fastmax, slow, slowmax (default: fast)')
ccd_group.add_argument('--seed',     metavar='INT',  type=int, help='Random number seed to bootstrap/reproduce results (default: random)')

hpc_group = parser.add_argument_group('PERFORMANCE')
hpc_group.add_argument('--data', metavar='PATH', type=str, help='Path RV file created by "rv-generator.py"')
hpc_group.add_argument('--dex',  metavar='NAME', type=int, help='Pyxel parallisation: index for running Pyxel on CPUs')
hpc_group.add_argument('--cpu',  metavar='INT',  type=str, help='PyEchelle parallisation: Number of CPUs for order-wise parallelisation')
hpc_group.add_argument('--cuda', action='store_true',      help='PyEchelle parallisation: Flag to use CUDA NVIDIA hardware for raytracing (makes cpu flag obsolete')
hpc_group.add_argument('--zip',  action='store_true',      help='Flag to zip output files.')

args = parser.parse_args()

# Create instance of class
m = marvelsim()
m.init_pyechelle()
pyxel_path = m.init_pyxel()

# Simulate a single frame

if args.frame is not None or args.debug is not None:

    # Set imgtype correctly
    if args.debug: args.frame = args.debug
    
    # Run bias or dark with numpy
    if args.frame == 'BBBBB' or args.frame == 'DDDDD':
        m.run_bias(args.frame)
        m.run_pyxel(args.frame)
        
    # Run everything else with PyEchelle and Pyxel
    else:
        m.run_pyechelle(args.frame)
        m.run_pyxel(args.frame)

# Simulate a set of calibrated spectra
        
elif args.calibs:
    # Create bias and darks fits with Numpy
    m.run_bias('BBBBB')
    m.run_bias('DDDDD')
    # Run pyechelle
    m.run_pyechelle('FFFFF')
    m.run_pyechelle('TTTTT')
    m.run_pyechelle('ETTTT')
    # Run Pyxel
    m.run_pyxel('BBBBB')
    m.run_pyxel('DDDDD')
    m.run_pyxel('FFFFF')
    m.run_pyxel('TTTTT')
    m.run_pyxel('ETTTT')

# Simulate science spectra
    
else:

    # Run pyechelle and Pyxel together
    if args.science:
        m.run_pyechelle('TSSSS')
        m.run_pyxel('TSSSS')

    # Run pyechelle alone with either CUDA or CPUs 
    elif args.cuda or args.cpu:
        m.run_pyechelle('ESSSS')

    # Run pyxel alone with CPUs
    elif args.dex:
        m.run_pyxel_cpu('ESSSS')
        #os.rmdir(pyxel_path)

    else:
        errorcode('error', 'Not valid setup for science mode! Use --science if attended')

# Final execution time
os.rmdir(pyxel_path)
toc = datetime.datetime.now()
print(f"\nMARVEL simulations took {toc - tic}")
