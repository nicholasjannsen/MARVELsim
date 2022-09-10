#!/usr/bin/env python3

"""
This code explores the python library RadVel to generate synthetic
Radial Velocity (RV) time series to be used as input for MARVELsim.
The RV input can be used to parallelise each MARVELsim stellar spectrum
using the High Performance Computer (HPC) facility VSC. Note that
a noise-less RV time series can be generated for multiple planets
if including corresponding amount of input values for each parameter.

User examples:
  $ python rv-generator.py  0 50 0 0 <path/to/output/file.txt>
  $ python rv-generator.py  10.0 30.0 0.0 150.0 <path/to/output/file.txt>

radvel documentation:
https://radvel.readthedocs.io/en/latest/tutorials/CustomModel-tutorial.html
"""

import os
import argparse
import radvel
import numpy as np
import matplotlib.pyplot as plt
from colorama import Fore, Style
from astropy import constants as c
from astropy import units as u
from utilities import errorcode

import matplotlib
matplotlib.rcParams.update({'font.size': 16})
plt.rcParams['text.usetex'] = True

#==============================================================#
#                           UTILITIES                          #
#==============================================================#

class RV(object):
    """
    This class has the purpose
    """
    
    # INITILIZE THE CLASS:
    def __init__(self):
        """
        Constructor of the class.
        """
        self.tdur = args.tdur


        
    
    def rv_model(self, t, t0, Ms, Mp, P, a, e, i, w):

        # Time of periastron
        tp = t0.to('d').value - P.to('d').value * (np.pi/2. - w.value)

        # True anomaly with radvel
        nu = radvel.orbit.true_anomaly(t, tp, P.to('d').value, e) * u.rad

        # RV signal as function of nu: Murray & Correia (2011) Eq. 61, 65 and 66
        # NOTE "astar" in the following is the reduced semimajor axies due to the common
        # center-of-mass and "K" is the relative RV semi-amplitude of the star.
        astar = Mp / (Mp + Ms) * a
        K     = 2.*np.pi/P * astar * np.sin(i)/np.sqrt(1. - np.power(e,2))
        RV    = K * (np.cos(nu + w) - e*np.cos(w))

        return RV, K





    def combine_model(self):
        
        # Prepare time points
        self.t  = np.linspace(0, self.tdur, self.tdur)
        self.tt = np.linspace(0, self.tdur, self.tdur*1000)

        # If constant star is requested
        if args.constant:
            self.RV = np.ones(len(self.t)) * args.constant

        # If RV light curve is requested
        else:

            # Single entry parameters
            self.Rs = (args.rs * u.R_sun).to('m')
            self.Ms = (args.ms * u.M_sun).to('kg')

            # Multi entry parameters
            self.t0 = (args.t0[0] * u.d).to('s')
            self.P  = (args.p[0] * u.d).to('s')
            self.e  = args.e[0] 
            self.i  = (args.i[0] * u.deg).to('rad')
            self.w  = (args.w[0] * u.deg).to('rad')
            self.Rp = (args.rp[0] * u.R_earth).to('m')
            self.Mp = (args.mp[0] * u.M_earth).to('kg')

            # Calculate semi-major axis (K3)
            self.a  = (c.G * self.P**2 * (self.Ms + self.Mp) / (4*np.pi**2))**(1/3.)

            # Count number of planets
            self.n_planets = len(self.t0)

            RV0 = []
            RV1 = []
            K   = []

            for n in range(self.n_planets):

                # Make models
                rv0, K0 = self.rv_model(self.t,  self.t0[n], self.Ms, self.Mp[n], self.P[n], self.a[n], self.e[n], self.i[n], self.w[n])
                rv1, _  = self.rv_model(self.tt, self.t0[n], self.Ms, self.Mp[n], self.P[n], self.a[n], self.e[n], self.i[n], self.w[n])

                RV0.append(rv0)
                RV1.append(rv1)
                K.append(f'{K0.value:.2f} m/s')


            # Combine models
            self.RV = np.sum(np.array(RV0), axis=0)
            self.rv = np.sum(np.array(RV1), axis=0)
            self.K  = ', '.join(K) 

            # Plot light curve if required
            self.plot_model()




    def save_model(self):

        # Select output format
        if args.constant:
            fmt = '%i'
        else:
            fmt = ['%i', '%0.6f', '%0.6f']
                    
        # We here save the output in a format that HPC worker can read
        if args.outputfile:
            index  = np.arange(1, len(self.RV)+1, 1)
            header = 'index, time, rv'
            np.savetxt(args.outputfile, np.transpose([index, self.t, self.RV]),
                       fmt=fmt, header=header, comments='', delimiter=',')
            # Save as feather instead
            #import pandas as pd
            #df = pd.DataFrame(np.transpose([self.t, self.RV]), columns=['time', 'rv0'])
            #df.to_feather(args.outputfile)



    def plot_model(self):

        # Plot and save model
        plt.figure(figsize=(10,5))

        # Prepare title
        Rs = self.Rs.to('R_sun').value
        Ms = self.Ms.to('M_sun').value

        lab_star = (r'\textbf{Star:}' +
                    r' $R_s$ = '+f'{Rs:.2f}'+r' $R_{\odot}$;' +
                    r' $M_s$ = '+f'{Ms:.2f}'+r' $M_{\odot}$')

        if self.n_planets == 1:

            # Revert parameters for plot
            Rp = self.Rp.to('R_earth').value[0]
            Mp = self.Mp.to('M_earth').value[0]
            t0 = self.t0.to('d').value[0]
            P  = self.P.to('d').value[0]
            a  = self.a.to('R_sun').value[0]
            i  = self.i.to('deg').value[0]
            w  = self.w.to('deg').value[0]
            e  = self.e[0]

            lab_planet = (r'\textbf{Planet:}' +
                          r' $R_p$ = '+f'{Rp:.1f}'+r' $R_{\oplus}$;' +
                          r' $M_p$ = '+f'{Mp:.1f}'+r' $M_{\oplus}$;' +
                          r' $t_0$ = '+f'{t0:.1f}'+r' days;' +
                          r' $P$ = '+f'{P:.1f}'+r' days;' +
                          r' $a$ = '+f'{a:.1f}'+r' $R_{\odot}$;' +
                          r' $i$ = '+f'{i:.1f}'+r'$^{\circ}$;' +
                          r' $w$ = '+f'{w:.1f}'+r'$^{\circ}$;' +
                          r' $e$ = '+f'{e:.1f}')

            plt.title(lab_star + '\n' + lab_planet, fontsize=14)
        else:
            plt.title(lab_star, fontsize=14)

        plt.plot(self.t,  self.RV, 'mo', alpha=0.5, label='Observations')
        plt.plot(self.tt, self.rv, 'k:', alpha=0.5, label=r'$K_{RV}$ = ' + self.K)
        plt.xlabel('Time [d]')
        plt.ylabel('Radial Velocity [m/s]')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.show()


#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-o', '--outputfile', metavar='PATH', type=str, help='Output file -> /path/to/filename.txt')

obs_group = parser.add_argument_group('OBSERVATION')
obs_group.add_argument('-tdur',     metavar='DAY', type=int, help='Duration of observing campaign [days]')
obs_group.add_argument('-constant', metavar='RV',  type=int, help='Flag to produce a constant light curve')

star_group = parser.add_argument_group('STAR')
star_group.add_argument('-rs', metavar='RSUN', type=float, help='Stellar radius [Rsun]')
star_group.add_argument('-ms', metavar='MSUN', type=float, help='Stellar mass   [Msun]')

planet_group = parser.add_argument_group('PLANET')
planet_group.add_argument('-t0', metavar='DAY',    action='append', type=float, nargs='*', help='Known ephemeris time(s) of observation [days]')
planet_group.add_argument('-p',  metavar='DAY',    action='append', type=float, nargs='*', help='Orbital period(s) of planet(s) [days]')
planet_group.add_argument('-i',  metavar='DEG',    action='append', type=float, nargs='*', help='Inclination of planets(s) orbits [deg]')
planet_group.add_argument('-w',  metavar='DEG',    action='append', type=float, nargs='*', help='Argument(s) of perioastron [deg]')
planet_group.add_argument('-e',  metavar='NUM',    action='append', type=float, nargs='*', help='Eccentrity of planets(s)')
planet_group.add_argument('-rp', metavar='REARTH', action='append', type=float, nargs='*', help='Planet(s) radius [Rearth]')
planet_group.add_argument('-mp', metavar='MEARTH', action='append', type=float, nargs='*', help='Planet(s) mass   [Mearth]')

args = parser.parse_args()

# Initialise instance of class

rv = RV()
rv.combine_model()
rv.save_model()


