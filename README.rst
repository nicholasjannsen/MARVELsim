README for simulating realistic spectra for MARVEL
==================================================

Files from Julian Stuermer (2021-11-25)
---------------------------------------

MARVEL_2021_11_22.hdf: MARVEL optical model incl. CCD efficiency
examples.md:           command line examples to simulate various spectra with pyechelle

Files from Jake Pember (2021-11-16)

MARVEL_2021_11_16.hdf: MARVEL optical configuration file for PyEchelle
MARVEL_2021_11_16.ZDA: Zemax file used to create the HDF file 
MARVEL_2021_11_16.zmx: Zemax file used to create the HDF file
MARVEL_2021_11_16.zar: Zemax archive file

transmission_data.csv: total estimated throughput on the detector as a function of wavelength. No QE accounted for.
                       numbers not final as this is highly dependent on the optical coatings, in particular for the collimator mirror (3 passes) 
                       and the Echelle grating. Supplier has not committed to a minimum efficiency.

transmission.py:       python code to read the transmission_data.csv file

Transmission:          folder containing transmission reports
|_ c??.txt:            transmission report of 12 configuration (echelle orders equally spaced in 'm' number). 
                       These break down the transmission for each of the 5 x 12 wavelengths into the loss at every surface.

MARVEL simulation examples
--------------------------

The commands shown below are supposed to be executed in the root directory of pyechelle.

Bias frame
..........

The bias level and std read noise value is measured from HERMES data. The values are bias = 2170 ADU and read-noise-std = 5.5 ADU, respectively. Since we need the count in electrons we multiply with the gain of 9.4 e/ADU which gives 20398 e and 52 e. 

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant -t 0.0 --bias 20398 --read_noise 52 -o marvel_bias.fits


Flat frame
..........

We need to investigate the flat exposure time

.. code-block:: shell
		
   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant Constant Constant Constant -t 10 --bias 20398 --read_noise 52 -o marvel_flat.fits

Etalon frame
............

We need to investigate the flat exposure time

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Etalon --etalon_d=6 -t 10 --bias 20398 --read_noise 52 -o marvel_flat.fits

Generate a science frame

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Phoenix Phoenix Phoenix Phoenix Etalon --etalon_d=6 --d_primary 0.8 --d_secondary 0.1 --phoenix_t_eff 5800 --phoenix_log_g 4.5 --phoenix_z 0.0 --phoenix_alpha 0.0 --phoenix_magnitude 10.0 -t 1200 -o output/marvel_science_G2V_10mag_1200s.fits

MARVELsim
---------


