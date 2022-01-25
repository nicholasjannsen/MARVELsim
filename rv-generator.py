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
from astropy import constants as c
from astropy import units as u
from utilities import errorcode

#==============================================================#
#                           UTILITIES                          #
#==============================================================#

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

parser.add_argument('-t0', action='append', type=float, nargs='*', help='Known ephemeris time(s) of observation [d]')
parser.add_argument('-p',  action='append', type=float, nargs='*', help='Orbital period(s) of planet(s) [d]')
parser.add_argument('-i',  action='append', type=float, nargs='*', help='Inclination of planets(s) orbits')
parser.add_argument('-e',  action='append', type=float, nargs='*', help='Eccentrity of planets(s)')
parser.add_argument('-w',  action='append', type=float, nargs='*', help='Argument(s) of perioastron [deg]')
parser.add_argument('-rs',  action='append', type=float, nargs='*', help='Stellar radius [Rsun]')
parser.add_argument('-ms',  action='append', type=float, nargs='*', help='Stellar mass   [Msun]')
parser.add_argument('-rp',  action='append', type=float, nargs='*', help='Planet(s) radius [Rsun]')
parser.add_argument('-mp',  action='append', type=float, nargs='*', help='Planet(s) mass   [Msun]')

parser.add_argument('-o',  type=str, help='Output directory and filename -> /path/to/filename.txt')

args = parser.parse_args()
t0 = args.t0[0]
P  = args.p[0]
i  = args.i[0]
e  = args.e[0]
w  = args.w[0]
Rs = args.rs[0]
Ms = args.ms[0]
Rp = args.rp[0]
Mp = args.mp[0]

n_planets = len(t0)

#print(Mp)
#exit()

#==============================================================#
#                         CUSTOM MODELS                        #
#==============================================================#

# Set defautls
if args.t0 is None: args.t0 = 2456300

# Calculate semi-major axis (K3)
p  = (P * u.d).to('s')
ms = (Ms * u.M_sun).to('kg')
mp = (Ms * u.M_earth).to('kg')
a  = ((c.G * p**2 * (ms + mp) / (4*np.pi**2))**(1/3.)).to('AU').value

# Initialize and setup radvel model 
params = radvel.Parameters(num_planets=n_planets)

# Setup model
xx = []
for n in range(n_planets):
    m = n + 1
    params[f'tc{m}']  = radvel.Parameter(value=t0[n])
    params[f'per{m}'] = radvel.Parameter(value=P[n])
    params[f'a{m}']   = radvel.Parameter(value=a[n])
    params[f'rp{m}']  = radvel.Parameter(value=Rp[n]) # Rs/Rp
    params[f'inc{m}'] = radvel.Parameter(value=i[n])
    params[f'e{m}']   = radvel.Parameter(value=e[n], vary=False, linear=False)
    params[f'w{m}']   = radvel.Parameter(value=w[n], vary=False, linear=False)
    params[f'k{m}']   = radvel.Parameter(value=30)

    x = 10 * n
    indices = {
        f'tc{m}': 0+x,
        f'per{m}': 1+x,
        f'rp{m}': 2+x,
        f'a{m}': 5+x,
        f'inc{m}': 4+x,
        f'e{m}': 5+x,
        f'w{m}': 6+x,
        f'k{m}': 7+x,
        f'dvdt': 8+x,
        f'curv': 9+x,
    }
    xx.append(indices)

if n_planets==1: indices = {**xx[0]}
if n_planets==2: indices = {**xx[0], **xx[1]}
if n_planets==3: indices = {**xx[0], **xx[1], **xx[2]}
if n_planets==4: indices = {**xx[0], **xx[1], **xx[2], **xx[3]}
if n_planets==5: indices = {**xx[0], **xx[1], **xx[2], **xx[3], **xx[4]}

#print(new); exit()

# indices = {
#     'tc1': 0,
#     'per1': 1,
#     'rp1': 2,
#     'a1': 3,
#     'inc1': 4,
#     'e1': 5,
#     'w1': 6,
#     'k1': 7,
#     'dvdt': 8,
#     'curv': 9,
#     'tc2': 10,
#     'per2': 11,
#     'rp2': 12,
#     'a2': 13,
#     'inc2': 14,
#     'e2': 15,
#     'w2': 16,
#     'k2': 17,
#     'dvdt': 18,
#     'curv': 19,

# }


    
# params = radvel.Parameters(num_planets=1)
# params['tc1']  = radvel.Parameter(value=t0[0])
# params['per1'] = radvel.Parameter(value=P[0])
# params['a1']   = radvel.Parameter(value=a[0])
# params['rp1']  = radvel.Parameter(value=Rp[0]) # Rs/Rp
# params['inc1'] = radvel.Parameter(value=i[0])
# params['e1']   = radvel.Parameter(value=e[0], vary=False, linear=False)
# params['w1']   = radvel.Parameter(value=w[0], vary=False, linear=False)
# params['k1']   = radvel.Parameter(value=30)


# Run RV model
mod_rv = radvel.GeneralRVModel(params, forward_model=rv_calc)
mod_rv.vector.indices = indices
mod_rv.vector.dict_to_vector()

# Plot and save model
t_rv = np.linspace(2360, 2470, 1000)
plt.figure(figsize=(12,6))
plt.plot(t_rv, mod_rv(t_rv), 'm-')
plt.xlabel('Time [d]')
plt.ylabel('Radial Velocity [m/s]')
plt.title('RV model')
plt.show()
