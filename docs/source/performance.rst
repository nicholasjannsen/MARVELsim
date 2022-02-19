.. _performance:

Performance
===========

In order to speed up the simulations in order obtain a set of calibrated data or a time series of stellar spectra with RV variations, we here demonstrate how to run MARVELsim on a High Performace Computing (HPC) facility. We note that MARVELsim is specifically designed to run with CUDA NVIDIA hardware which typically is available for GPU nodes on most computing clusters. Since PyEchelle are the bottleneck for running large simulations using CUDA or running the simulations on multiple CPUs has a dramatic decrease of the computational time. We specifically refer to `PyEchelle's documentation on performance <https://stuermer.gitlab.io/pyechelle/benchmark.html>`_ for more information.

We will use the Vlaams Supercomputing Centre (VSC) as example on how to create a job script that runs with GPUs. We likewise provide two job script examples in the ``MARVELsim/hpc`` folder for your convenience. The first (``run_calibs.pbs``) shows a typical job script for running a full set of calibrated data:

.. code-block:: shell

   #!/bin/bash

   #PBS -N output
   #PBS -A <account_name>
   #PBS -l nodes=1:ppn=36:gpus=4:skylake
   #PBS -l partition=gpu
   #PBS -l pmem=2gb
   #PBS -l walltime=03:00:00

   cd $PBS_O_WORKDIR

   PYTHONPATH=$VSC_DATA/MARVELsim/marvelsim/bin/python
   export PYTHONPATH
   SIMDIR=$VSC_DATA/MARVELsim
   export SIMDIR

   # Activate environment 
   source marvelsim/bin/activate

   # Run MARVELsim
   cd $SIMDIR
   python simulator-marvel.py --calibs --cuda --zip -o $SIMDIR/output

Illustrated here we request a single node with 4 GPUs using each using 9 CPU claves (hence 36 in total) to execute the job. We request 2 GB of memory RAM to be on the safe side since a single 10,560 x 10,560 pixel full frame image occupy 851 Mb. Notice that it is possible to compress each image on the fly by enabling the flag ``zip`` as done in this example. Typical deflation rates per image are around 80%, hence, it is highly recommended to invoke this flag for faster data transfer after end job. For the job script show above the total run time (walltime) was 2 hours and 40 minutes.

The following example (``run_science.pbs``) shows a job script for running 300 stellar spectra using a generated RV time series called ``rv_data.txt``:

.. code-block:: shell

   #!/bin/bash

   #PBS -N output
   #PBS -A <account/project>
   #PBS -l nodes=1:ppn=36:gpus=4:skylake
   #PBS -l partition=gpu
   #PBS -l pmem=2gb
   #PBS -l walltime=40:00:00

   cd $PBS_O_WORKDIR

   PYTHONPATH=$VSC_DATA/MARVELsim/marvelsim/bin/python
   export PYTHONPATH
   SIMDIR=$VSC_DATA/MARVELsim
   export SIMDIR

   # Activate environment 
   source marvelsim/bin/activate

   # Run MARVELsim
   cd $SIMDIR
   python simulator-marvel.py --time 900 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --data rv_data.txt --cuda --zip -o $SIMDIR/output

Compared to the previous job script we here use the same computational resources, however, with the exception of increasing the walltime. Notice that adding more nodes will not speed up the computations, however, some cluster do provide more GPUs which will decrease the run time. We recommend to debug and test the computational resources needed for your jobs adding ``#PBS -l qos=debugging`` to the PSB details in the scripts shown above.  
