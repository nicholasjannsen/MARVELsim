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

def rv_model(t, t0, Ms, Mp, P, a, e, i, w):

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


#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-o', metavar='OUTFILE', type=str, help='Output file -> /path/to/filename.txt')

obs_group = parser.add_argument_group('OBSERVATION')
obs_group.add_argument('-tdur', metavar='DAY', type=int, help='Duration of observing campaign [days]')

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

# Single entry parameters
tdur = args.tdur
Rs = (args.rs * u.R_sun).to('m')
Ms = (args.ms * u.M_sun).to('kg')

# Multi entry parameters
t0 = (args.t0[0] * u.d).to('s')
P  = (args.p[0] * u.d).to('s')
e  = args.e[0] 
i  = (args.i[0] * u.deg).to('rad')
w  = (args.w[0] * u.deg).to('rad')
Rp = (args.rp[0] * u.R_earth).to('m')
Mp = (args.mp[0] * u.M_earth).to('kg')

# Calculate semi-major axis (K3)
a  = (c.G * P**2 * (Ms + Mp) / (4*np.pi**2))**(1/3.)

# Count number of planets
n_planets = len(t0)

#==============================================================#
#                         CUSTOM MODELS                        #
#==============================================================#

# Prepare time points
t  = np.linspace(0, tdur, tdur)
tt = np.linspace(0, tdur, tdur*1000)

RV0 = []
RV1 = []
K   = []

for n in range(n_planets):

    # Make models
    rv0, K0 = rv_model(t,  t0[n], Ms, Mp[n], P[n], a[n], e[n], i[n], w[n])
    rv1, _  = rv_model(tt, t0[n], Ms, Mp[n], P[n], a[n], e[n], i[n], w[n])

    RV0.append(rv0)
    RV1.append(rv1)
    K.append(f'{K0.value:.2f} m/s')


# Combine models
RV = np.sum(np.array(RV0), axis=0)
rv = np.sum(np.array(RV1), axis=0)
K  = ', '.join(K) 

#==============================================================#
#                           MAKE PLOT                          #
#==============================================================#

# Plot and save model
plt.figure(figsize=(10,5))

# Prepare title
Rs = Rs.to('R_sun').value
Ms = Ms.to('M_sun').value

lab_star   = (r'\textbf{Star:}' +
              r' $R_s$ = '+f'{Rs:.2f}'+r' $R_{\odot}$;' +
              r' $M_s$ = '+f'{Ms:.2f}'+r' $M_{\odot}$')

if n_planets == 1:
    
    # Revert parameters for plot
    Rp = Rp.to('R_earth').value[0]
    Mp = Mp.to('M_earth').value[0]
    t0 = t0.to('d').value[0]
    P  = P.to('d').value[0]
    a  = a.to('R_sun').value[0]
    i  = i.to('deg').value[0]
    w  = w.to('deg').value[0]
    e  = e[0]

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
    
plt.plot(t,  RV, 'mo', alpha=0.5, label='Observations')
plt.plot(tt, rv, 'k:', alpha=0.5, label=r'$K_{RV}$ = ' + K)
plt.xlabel('Time [d]')
plt.ylabel('Radial Velocity [m/s]')
plt.legend(loc='best')
plt.tight_layout()
plt.show()





# FUTRUE WORK

# # Initialize and setup radvel model 
# params = radvel.Parameters(num_planets=n_planets)

# # Setup model
# xx = []
# for n in range(n_planets):
#     m = n + 1
#     params[f'tc{m}']  = radvel.Parameter(value=t0[n])
#     params[f'per{m}'] = radvel.Parameter(value=P[n])
#     params[f'a{m}']   = radvel.Parameter(value=a[n])
#     params[f'rp{m}']  = radvel.Parameter(value=Rp[n]) # Rs/Rp
#     params[f'inc{m}'] = radvel.Parameter(value=i[n]*np.pi/180)
#     params[f'e{m}']   = radvel.Parameter(value=e[n], vary=False, linear=False)
#     params[f'w{m}']   = radvel.Parameter(value=w[n]*np.pi/180, vary=False, linear=False)
#     params[f'k{m}']   = radvel.Parameter(value=30, vary=False)

#     x = 10 * n
#     indices = {
#         f'tc{m}': 0+x,
#         f'per{m}': 1+x,
#         f'rp{m}': 2+x,
#         f'a{m}': 5+x,
#         f'inc{m}': 4+x,
#         f'e{m}': 5+x,
#         f'w{m}': 6+x,
#         f'k{m}': 7+x,
#         f'dvdt': 8+x,
#         f'curv': 9+x,
#     }
#     xx.append(indices)

# if n_planets==1: indices = {**xx[0]}
# if n_planets==2: indices = {**xx[0], **xx[1]}
# if n_planets==3: indices = {**xx[0], **xx[1], **xx[2]}
# if n_planets==4: indices = {**xx[0], **xx[1], **xx[2], **xx[3]}
# if n_planets==5: indices = {**xx[0], **xx[1], **xx[2], **xx[3], **xx[4]}

# #print(params); exit()

# # Run RV model
# mod_rv = radvel.GeneralRVModel(params, forward_model=rv_calc)
# mod_rv.vector.indices = indices
# mod_rv.vector.dict_to_vector()
