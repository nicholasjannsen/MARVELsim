.. _extra examples:

Extra examples
==============

This section serves as an extra overview of the functionality of PyEchelle and Pyxel, respectively. 

PyEchelle examples
------------------

The following examples shows the ease of using PyEchelle from bash. Notice that if you run ``pyechelle`` from a folder different that the base of MARVELsim, a new ``model`` folder will created in order to download the spectrograph model. Here we show a few example while using model of the MARVEL spectrograph:

**Bias images:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant -t 0.0 --bias 2000 --read_noise 5 -o marvel_bias.fits


**Spectral flats:**

.. code-block:: shell
		
   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant Constant Constant Constant -t 5 --bias 20398 --read_noise 52 -o marvel_flat.fits

**Etalon spectra:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Etalon --etalon_d=6 -t 10 --bias 20398 --read_noise 52 -o marvel_flat.fits

**Stellar spectra:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Phoenix Phoenix Phoenix Phoenix Etalon --etalon_d=6 --d_primary 0.8 --d_secondary 0.1 --phoenix_t_eff 5800 --phoenix_log_g 4.5 --phoenix_z 0.0 --phoenix_alpha 0.0 --phoenix_magnitude 10.0 -t 1200 -o output/marvel_science_G2V_10mag_1200s.fits

.. warning::

   Note that the spectral flat simulates a large amount of photons and thus is the most time consuming. We strongly suggest to use a supported NVIDIA driver for CUDA or split the simulations on several CPUs for parallel computing. Please have a look at the :ref:`Performance section <performance>`.
   

Pyxel examples
--------------



