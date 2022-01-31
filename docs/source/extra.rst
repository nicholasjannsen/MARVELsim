Extra examples
==============

PyEchelle examples
------------------

The commands shown below are supposed to be executed in the root directory of pyechelle.

**Bias images:** The bias level and std read noise value is measured from HERMES data. The values are bias = 2170 ADU and read-noise-std = 5.5 ADU, respectively. Since we need the count in electrons we multiply with the gain of 9.4 e/ADU which gives 20398 e and 52 e. 

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant -t 0.0 --bias 20398 --read_noise 52 -o marvel_bias.fits


**Spectral flats:**

.. code-block:: shell
		
   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant Constant Constant Constant -t 10 --bias 20398 --read_noise 52 -o marvel_flat.fits

**Etalon spectra:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Etalon --etalon_d=6 -t 10 --bias 20398 --read_noise 52 -o marvel_flat.fits

**Stellar spectra:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Phoenix Phoenix Phoenix Phoenix Etalon --etalon_d=6 --d_primary 0.8 --d_secondary 0.1 --phoenix_t_eff 5800 --phoenix_log_g 4.5 --phoenix_z 0.0 --phoenix_alpha 0.0 --phoenix_magnitude 10.0 -t 1200 -o output/marvel_science_G2V_10mag_1200s.fits


Pyxel examples
--------------



