Installation
============

In this tutorial we cover how to setup an Python environment that allows you quickly get started to run simulations of MARVEL spectra.

.. warning::

   At the moment this software is only supported for Python 3.8.

   
1. Downlod from source
----------------------

Move to a desired directory for which you want to download the software and simply:

.. code-block:: shell

   git clone https://github.com/nicholasjannsen/MARVELsim.git


2. Create Conda environment
---------------------------
   
We strongly recommend to install Poetry through a Anaconda environment. This is especially handy when installing on a computing cluster since these typically only have limited versions of Python installed by default. While using Anaconda all Python versions can be installed and thus an exact freeze of the poetry installation can be made. Thus first install Anaconda or miniconda and the create a new virtuel environment called ``marvelsim``: 

.. code-block:: shell
		
   conda create -n marvelsim python=3.8 numpy scipy matplotlib

Now activate your new conda environment:

.. code-block:: shell

   conda activate marvelsim

   
3. Install with Poetry
----------------------

Since the most commen use case for MARVELsim is to run with HPC, we use `Poetry <https://python-poetry.org/>`_ to manage and install out Python libraries. First `install Poetry <https://python-poetry.org/docs/master/>`_ from the **master** branch 

.. code-block:: shell

   curl -sSL https://install.python-poetry.org | python -
   
Verify that poetry was installed successfully by typing

``poetry --version``

We commanded to include the following path to your ``~/.bashrc`` file:

.. code-block:: shell

   export PATH="$HOME/.poetry/bin:$PATH"

Deactivate your conda environment and install with poetry:

.. code-block:: shell

   conda deactivate
   poetry install

   
3. Install with Poetry on HPC
-----------------------------

Here we show how to install Poetry on the VSC. When installing on a computing cluster you should typically want to avoid installing software directly to your ``$HOME`` space since this has a very limited amount of available space. Must clusters recommend to install software in your ``$DATA`` directory. Hence:

.. code-block:: shell

   curl -sSL https://install.python-poetry.org | POETRY_HOME=$VSC_DATA/poetry python -
   
Verify that poetry was installed successfully by typing ``poetry --version``. In order for Poetry to available from any compute node you need to include the following path to your ``~/.bashrc`` file:

.. code-block:: shell

   POETRY=$VSC_DATA/poetry/bin
   export POETRY

Next change the installation location of the virtuel poetry environment to:
   
.. code-block:: shell

   poetry config virtualenvs.path $VSC_DATA/poetry/virtualenvs

Finally deactiavte your conda environment and install MARVELsim from the base directory using:

.. code-block:: shell

   conda deactivate
   poetry install

.. warning::

   It is important that the conda environment is deactivated before running ``poetry install``. If not done all packages will be installed within your conda ``bin`` folder and poetry cannot find your packges when using either ``poetry shell`` to spawn a shell environment or ``poetry run python marvelsim.py`` for running a script directly.
  

Install with a python3-venv
---------------------------
   
Another method is use a control all dependencies with a python3-venv environment. First install the python3-venv package

.. code-block:: shell

   sudo apt install python3-venv

Next create a virtual environment called ``marvelsim`` within the base of the cloned repositry

.. code-block:: shell

   python -m venv marvelsim

To activate and deactivate the environment simply use

.. code-block:: shell
		
   source marvelsim/bin/activate
   deactivate

   
Extra tools
-----------

Before starting investigating your output fits files we recomment to install `dfits <https://www.eso.org/sci/software/eclipse/eug/eug/node8.html>`_ which is an nice tool to inspect fits headers (e.g. ``dfits <filename>.fits``). On Linux install this packge with:

.. code-block:: shell

   sudo apt-get install qfits-tools

In addition, the astronomy software `ds9 <https://sites.google.com/cfa.harvard.edu/saoimageds9>`_ is an indispensable tool to quickly view your fits images (e.g. ``ds9 <filename>.fits``). Install this software with:

.. code-block:: shell

   sudo apt install saods9
