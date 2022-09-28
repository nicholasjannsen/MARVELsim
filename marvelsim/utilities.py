#!/usr/bin/env python3

"""
Placeholder for all the small utilities in this repo.
"""

import pkg_resources
from astropy.io import fits
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



def add_fitsheader(filename, obsmode, exptime,
                   readmode, speed, gain_ac, bias,
                   readsen, gain_amp, adc_range):
    # Load and open file
    hdul = fits.open(filename)
    hdr = hdul[0].header
    # Get software versions
    marvelsim_version = pkg_resources.get_distribution('MARVELsim').version
    pyechelle_version = pkg_resources.get_distribution('pyechelle').version
    # Add headers (pyxel version is added automatically)
    hdr.append(('PYECHE_V', f'{pyechelle_version}', 'PyEchelle version'), end=True)
    hdr.append(('MARVEL_V', f'{marvelsim_version}', 'MARVELsim version'), end=True)
    hdr.append(('FILENAME', f'{filename}', 'Origina'), end=True)
    hdr.append(('OBSMODE',  f'{obsmode}', 'Observing mode'), end=True)
    hdr.append(('IMAGETYP', f'{obsmode}', 'Image type'), end=True)
    hdr.append(('EXPTYPE',  f'{obsmode}', 'Exposure type'), end=True)
    #------------------------
    hdr.append(('ORIGIN', 'Instituut voor Sterrenkunde, KU Leuven', 'Institution'), end=True)
    hdr.append(('OBSERVAT', 'LaPalma', 'Observatory name'), end=True)
    hdr.append(('TELESCOP', 'MARVEL', 'Telescope name'), end=True)
    hdr.append(('OBSGEO-X', '5327306.5552', 'Cartesian X [meters] GRS80'), end=True)
    hdr.append(('OBSGEO-Y', '-1718448.6952', 'Cartesian Y [meters] GRS80'), end=True)
    hdr.append(('OBSGEO-Z', '3051947.7715', 'Cartesian Z [meters] GRS80'), end=True)
    hdr.append(('OBSERVER', 'Nicholas Jannsen', 'Observer'), end=True)
    hdr.append(('PROG_ID', '0', 'Programme ID'), end=True)
    hdr.append(('INSTRUME', 'MARVEL', 'Instrument'), end=True)
    hdr.append(('FIBMODE', 'High-Res', 'Fibre mode'), end=True)
    hdr.append(('CREATOR', 'revision_20201117', 'Version of data acquisition system'), end=True)
    hdr.append(('HDRVERS', '20151026', 'Version of FITS header'), end=True)
    hdr.append(('UNSEQ', '995604', 'Unique sequence number'), end=True)
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
    hdr.append(('PCIFILE', 'None', 'PCI card setup file'), end=True)
    hdr.append(('TIMFILE', '/home/mocs/mocs/config/mocs/marvel/tim-MARVEL-20130205-sp_idle.lod', ''), end=True)
    hdr.append(('UTILFILE', 'None', 'Utility board setup file'), end=True)
    #------------------------
    hdr.append(('READMODE', f'{readmode}', 'Detector readout mode'), end=True)
    hdr.append(('DETSPEED', f'{speed}', '[kHz] Controller readout speed'), end=True)
    hdr.append(('DETGAIN',  f'{gain_ac}', '[e-/ADU] Detector gain'), end=True)
    hdr.append(('DETBIAS',  f'{bias}', '[ADU] Expected bias level'), end=True)
    hdr.append(('READSEN',  f'{readsen}', '[V/e-] Charge readout sensitivity'), end=True)
    hdr.append(('GAINAMP',  f'{gain_amp}', '[V/V] Pre-amplifier gain'), end=True)
    hdr.append(('ADCRANGE', f'{adc_range}', '[V] Output DC level'), end=True)
    #------------------------
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
    hdr.append(('TEMPH047', '16.996', '[C] MARVEL temperature air.camera'), end=True)
    hdr.append(('TEMPH048', '17.851', '[C] MARVEL temperature table.center'), end=True)
    hdr.append(('TEMPH050', '17.920', '[C] MARVEL temperature grating.mount.top'), end=True)
    hdr.append(('TEMPH051', '17.819', '[C] MARVEL temperature fiberexit.mount'), end=True)
    hdr.append(('TEMPH061', '17.859', '[C] MARVEL temperature maincoll.glass.top'), end=True)
    hdr.append(('TEMPH052', '18.030', '[C] MARVEL temperature maincoll.mount.top'), end=True)
    hdr.append(('TEMPH053', '17.977', '[C] MARVEL temperature maincoll.mount.bot'), end=True)
    hdr.append(('TEMPH054', '17.670', '[C] MARVEL temperature camera.top.center'), end=True)
    hdr.append(('TEMPH055', '16.483', '[C] MARVEL temperature cryostat.front'), end=True)
    hdr.append(('TEMPH056', '16.025', '[C] MARVEL temperature cryostat.rear'), end=True)
    hdr.append(('TEMPH057', '17.859', '[C] MARVEL temperature grating.glass.center'), end=True)
    hdr.append(('TEMPH058', '17.966', '[C] MARVEL temperature maincoll.glass.top'), end=True)
    hdr.append(('TEMPH059', '17.979', '[C] MARVEL temperature maincoll.glass.bot'), end=True)
    hdr.append(('BSCALE', '1', ''), end=True)
    hdr.append(('BZERO', '32768', ''), end=True)
    # Write new file with header
    fits.writeto(filename, hdul[0].data, hdr, overwrite=True)
