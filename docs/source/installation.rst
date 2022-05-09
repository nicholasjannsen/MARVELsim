Installation
============

In this tutorial we cover how to setup an Python environment that allows you quickly get started to run simulations of MARVEL spectra.

.. important::

   We recommend to follow the installation instructions below and using `Poetry <https://python-poetry.org/>`_ to install MARVELsim. While using Python 3.8 and Poetry you will get a direct freeze of the software which successfully will build across any platform.

   
1. Downlod from source
----------------------

Move to a desired directory for which you want to download the software and simply:

.. code-block:: shell

   git clone https://github.com/nicholasjannsen/MARVELsim.git


2. Create Conda environment
---------------------------
   
We strongly recommend to install Poetry through a `Conda <https://docs.conda.io/en/latest/>`_ environment. This is especially handy when installing on a computing cluster since these typically only have limited versions of Python installed by default. While using Conda all Python versions can be installed and thus an exact freeze of the poetry installation can be made. Thus first install Conda (or *Miniconda*) and the create a new environment called ``marvelsim``: 

.. code-block:: shell
		
   conda create -n marvelsim python=3.8

Now activate your new conda environment:

.. code-block:: shell

   conda activate marvelsim

   
3. Local (Poetry) installation
------------------------------

Since the most commen use case for MARVELsim is to run with HPC, it is recommended to use `Poetry <https://python-poetry.org/>`_ to manage and install Python libraries. First `install Poetry <https://python-poetry.org/docs/master/>`_ from the **master** branch 

.. code-block:: shell

   pip install --user poetry
   
Verify that poetry was installed successfully by typing ``poetry --version``. Next include the following path to your ``~/.bashrc`` file:

.. code-block:: shell

   POETRY=$HOME/.poetry/bin
   export POETRY

Deactivate your conda environment (or any other active Python environment) and install with your new poetry environment:

.. code-block:: shell

   conda deactivate
   poetry install

   
3. Cluster (Poetry) installation
--------------------------------

Here we show how to install Poetry on the VSC, however, the general workflow can most likely be followed for installing MARVELsim on any cluster. On a computing cluster you typically want to avoid installing software directly to your ``$HOME`` workspace since this has a very limited amount of available storage memory. Thus, for clusters we recommend to install software in your ``$DATA`` directory. Hence, we need to tell where to install Poetry:

.. code-block:: shell

   conda activate marvelsim
   pip install --user poetry | POETRY_HOME=$VSC_DATA/poetry python -
   
Verify that poetry was installed successfully by typing ``poetry --version`` and verify that installation location with ``which poetry``. Next change the installation location of the virtuel poetry eironment to:
   
.. code-block:: shell

   poetry config virtualenvs.path $VSC_DATA/poetry/virtualenvs

In order for Poetry to be available from any compute node, you need to include the following path to your ``~/.bashrc`` file:
   
.. code-block:: shell

   POETRY=$VSC_DATA/poetry/bin
   export POETRY

Finally deactiavte your Conda environment and install MARVELsim from the base directory using:

.. code-block:: shell

   conda deactivate
   poetry install
   
   
Extra tools
-----------

Before starting investigating your output fits files we recomment to install `dfits <https://www.eso.org/sci/software/eclipse/eug/eug/node8.html>`_ which is an nice tool to inspect fits headers (e.g. ``dfits <filename>.fits``). On Linux install this packge with:

.. code-block:: shell

   sudo apt-get install qfits-tools

In addition, the astronomy software `ds9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_ is an indispensable tool to quickly view your fits images (e.g. ``ds9 <filename>.fits``). Install this software with:

.. code-block:: shell

   sudo apt install saods9
