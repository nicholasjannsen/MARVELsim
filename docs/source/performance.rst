.. _performance:

Performance
===========

In order to speed up the simulations in order obtain a set of calibrated data or a time series of stellar spectra with RV variations, we here demonstrate how to run MARVELsim on a High Performace Computing (HPC) facility. We note that MARVELsim is specifically designed to run with CUDA NVIDIA hardware which typically is available for GPU nodes on most computing clusters. Since PyEchelle are the bottleneck for running large simulations using CUDA or running the simulations on multiple CPUs has a dramatic decrease of the computational time. We specifically refer to `PyEchelle's documentation on performance <https://stuermer.gitlab.io/pyechelle/benchmark.html>`_ for more information.

We will use the Vlaams Supercomputing Centre (VSC)

You encounter that the software cannot be installed on a HPC cluster due to a froozen version of Python. PyEchelle can only be installed through Python >=3.8. Thus we here show how to use Anaconda to install all packages. Some packages cannot be installed directly through Anaconda and thus needs to be installed with pip your newly created conda environment.

First let's start creating the conda environment:

.. code-block:: shell

   
