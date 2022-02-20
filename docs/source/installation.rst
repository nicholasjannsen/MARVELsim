Installation
============

In this tutorial we cover how to setup an Python environment that allows you quickly get started to run simulations of MARVEL spectra.

.. warning::

   At the moment this software is only supported for Python 3.8.

1. Downlod from source
-------------------

Move to a desired directory for which you want to download the software and simply:

.. code-block:: shell

   git clone https://github.com/nicholasjannsen/MARVELsim.git


2. Create Python environment
-------------------------

The easiest way version control both PyEchelle, Pyxel, and their dependencies is to create a python environment. First install the python3-venv package

.. code-block:: shell

   sudo apt install python3-venv

Next create a virtual environment called ``marvelsim`` within the base of the cloned repositry

.. code-block:: shell

   python -m venv marvelsim

To activate and deactivate the environment simply use

.. code-block:: shell
		
   source marvelsim/bin/activate
   deactivate


3. Install software
----------------

At the root of the MARVELsim repository a installation script called ``install.sh`` are provided for the installation. This scripts takes care of activating your virtual environment before installing the necessary packages using ``pip``. Now install all necessary Python libraries by simply commanding:

.. code-block:: shell

   ./install.sh


.. warning::

   The install script automatically installs Pyxel version 11.5, however, a newer version not yet well document and not supported by MARVELsim are available. Thus, double check that your have the correct version of pyxel by ``pyxel --version``.


Extra tools
-----------

Before starting investigating your output fits files we recomment to install `dfits <https://www.eso.org/sci/software/eclipse/eug/eug/node8.html>`_ which is an nice tool to inspect fits headers (e.g. ``dfits <filename>.fits``). On Linux install this packge with:

.. code-block:: shell

   sudo apt-get install qfits-tools

In addition, the astronomy software `ds9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_ is an indispensable tool to quickly view your fits images (e.g. ``ds9 <filename>.fits``). Install this software with:

.. code-block:: shell

   sudo apt install saods9
