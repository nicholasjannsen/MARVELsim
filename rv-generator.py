#!/usr/bin/env python3

"""
This code explores the python library "radvel" to generate synthetic
Radial Velocity (RV) time series to be used as input for MARVELsim.
The RV input can be used to parallelise each MARVELsim stellar spectrum
using the High Performance Computer (HPC) facility VSC.

User examples:
  $ python rv-generator.py  0 50 0 0 output/rv_test0.txt
  $ python rv-generator.py  10.0 30.0 0.0 150.0 output/rv_test0.txt

radvel documentation:
https://radvel.readthedocs.io/en/latest/tutorials/CustomModel-tutorial.html
"""

import os
import argparse
import radvel
#import pymc3 as pm
#import exoplanet as xo
import numpy as np
import matplotlib.pyplot as plt
from colorama import Fore, Style

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


def rv_calc(t, params, vector):

    per = vector.vector[1][0]
    tp = radvel.orbit.timetrans_to_timeperi(tc=vector.vector[0][0], per=vector.vector[1][0],
                                            ecc=vector.vector[5][0], omega=vector.vector[6][0])
    e = vector.vector[5][0]
    w = vector.vector[6][0]
    k = vector.vector[7][0]
    orbel_synth = np.array([per, tp, e, w, k])
    vel = radvel.kepler.rv_drive(t, orbel_synth)

    return vel


#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-t0', action='append', type=float, nargs='*',
                    help='Known ephemeris time(s) of observation [d]')
parser.add_argument('-p', action='append', type=float, nargs='*',
                    help='Orbital period(s) of planet(s) [d]')
parser.add_argument('-e', action='append', type=float, nargs='*',
                    help='Eccentrity of planets(s)')
parser.add_argument('-w', action='append', type=float, nargs='*',
                    help='Argument(s) of perioastron [deg]')
parser.add_argument('-o', type=str, help='Output directory and filename -> /path/to/filename.txt')

args = parser.parse_args()

#==============================================================#
#                         CUSTOM MODELS                        #
#==============================================================#

# Initialize and setup radvel model 

n_planets = len(args.t0[0])

params = radvel.Parameters(num_planets=1)
# if no_planets == 1:
#     P 
# for

# Setup model

for n in range(n_planets):

    params['tc'] = radvel.Parameter(value=2456300)
    params['per'] = radvel.Parameter(value=200)
    params['a'] = radvel.Parameter(value=10)
    params['rp'] = radvel.Parameter(value=0.08)
    params['inc'] = radvel.Parameter(value=90)
    params['e'] = radvel.Parameter(value=0.0, vary=False, linear=False)
    params['w'] = radvel.Parameter(value=0.0, vary=False, linear=False)
    params['k'] = radvel.Parameter(value=30)

    # params['tc']  = radvel.Parameter(value=args.t0[0][n])
    # params['per'] = radvel.Parameter(value=args.p[0][n])
    # params['e']   = radvel.Parameter(value=args.e[0][n], vary=False, linear=False)
    # params['w']   = radvel.Parameter(value=args.w[0][n], vary=False, linear=False)
    # params['a']   = radvel.Parameter(value=10)
    # params['rp']  = radvel.Parameter(value=0.08)
    # params['inc'] = radvel.Parameter(value=90)
    # params['k']   = radvel.Parameter(value=30)
    #params['jit_rv'] = radvel.Parameter(value=1.0)
    #params['gamma_rv'] = radvel.Parameter(value=0.0)

indices = {
    'tc': 0,
    'per': 1,
    'rp': 2,
    'a': 5,
    'inc': 4,
    'e': 5,
    'w': 6,
    'k': 7,
    'dvdt': 8,
    'curv': 9,
    #'jit_rv': 12,
    #'gamma_rv': 13
}


# for n in range(n_planets):
#     params['tc']  = radvel.Parameter(value=args.t0[0][n])
#     params['per'] = radvel.Parameter(value=args.p[0][n])
#     params['e']   = radvel.Parameter(value=args.e[0][n], vary=False, linear=False)
#     params['w']   = radvel.Parameter(value=args.w[0][n], vary=False, linear=False)
#     params['jit_trans'] = radvel.Parameter(value=0.01)
#     params['gamma_trans'] = radvel.Parameter(value=0, vary=False) #Unless you construct your own lik

# indices = {
#     'tc': 0,
#     'per': 1,
#     'e': 5,
#     'w': 6,
#     'dvdt': 8,
#     'curv': 9,
#     'jit_trans': 10,
#     'gamma_trans':11,
#     'jit_rv': 12,
#     'gamma_rv': 13
# }


# Run RV model

mod_rv = radvel.GeneralRVModel(params, forward_model=rv_calc)
mod_rv.vector.indices = indices
mod_rv.vector.dict_to_vector()

# Plot and save model

t_rv = np.linspace(2456200, 2457140, 1000)
plt.figure(figsize=(12,6))
plt.plot(t_rv, mod_rv(t_rv), 'm-')
plt.xlabel('Time [d]')
plt.ylabel('Radial Velocity [m/s]')
plt.title('RV model')
plt.show()
