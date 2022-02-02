.. _extra examples:

Extra examples
==============

This section serves as an extra overview of the functionality of PyEchelle and Pyxel, respectively. 

PyEchelle examples
------------------

The following examples shows the ease of using PyEchelle from bash. Notice that if you run ``pyechelle`` from a folder different that the base of MARVELsim, a new ``model`` folder will created in order to download the spectrograph model. Here we show a few example while using model of the MARVEL spectrograph and we here include a bias level of 10,000 e- including a root-mean-square read noise of 25 e- (which are equivalent ot 2000 ADU and 5 ADU, respectively, with a CCD gain of 5 e-/ADU):

**Bias images:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant -t 0.0 --bias 10000 --read_noise 25 -o marvel_bias.fits

**Spectral flats:**

.. code-block:: shell
		
   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Constant Constant Constant Constant -t 5 --bias 10000 --read_noise 25 -o marvel_flat.fits

**ThAr or ThNe Arc:**

.. code-block:: shell
		
   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources ThAr ThAr ThAr ThAr -t 5 --bias 10000 --read_noise 25 -o marvel_thar.fits
   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources ThNe ThNe ThNe ThNe -t 5 --bias 10000 --read_noise 25 -o marvel_thne.fits
   
**Etalon spectra:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Etalon --etalon_d=6 -t 10 --bias 10000 --read_noise 25 -o marvel_flat.fits

**Stellar spectra:**

.. code-block:: shell

   pyechelle --fiber 1-5 -s MARVEL_2021_11_22 --sources Phoenix Phoenix Phoenix Phoenix Etalon --etalon_d=6 --d_primary 0.8 --d_secondary 0.1 --phoenix_t_eff 5800 --phoenix_log_g 4.5 --phoenix_z 0.0 --phoenix_alpha 0.0 --phoenix_magnitude 10.0 -t 1200 -o marvel_science.fits

.. warning::

   Note that the spectral flat simulates a large amount of photons and thus is the most time consuming. We strongly suggest to use a supported NVIDIA driver for CUDA or split the simulations on several CPUs for parallel computing. Please have a look at the :ref:`Performance section <performance>`.
   

Pyxel examples
--------------

Pyxel has a quite ellaborate documentation page, but to make it easier to test what happens if including different modules and CCD effect to the MARVEL spectra, the folder called ``test_pyxel`` serves as a small testbed. Wihtin there is a very simply Python script called ``test.py`` which includes the following code:

.. code-block:: python

   #!/usr/bin/env python3

   import pyxel

   # Load specific inputfile for MARVEL
   config = pyxel.load("test.yaml")

   # Setup configurations
   exposure = config.exposure        # class Observation
   detector = config.ccd_detector    # class CCD
   pipeline = config.pipeline        # class DetectionPipeline

   # Run pyxel in exposure mode
   pipeline.charge_generation.load_charge.arguments.filename   = 'science.fits'
   pipeline.charge_generation.load_charge.arguments.time_scale = 5
   pipeline.charge_generation.tars.arguments.seed              = np.random.randint(1e9)

   results = pyxel.exposure_mode(exposure=exposure, detector=detector, pipeline=pipeline)

Seen above scirpt simply initializes Pyxel and reads the ``test.yaml`` inputfile configuration file. We here here show an example of running Pyxel in exposure mode and changing a few parameters, but this can also easily be done directly from within the configuration yaml file. To lower the computation time, running this test script will return a Pyxel fits images of 300x300 pixels.  
