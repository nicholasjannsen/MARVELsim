.. _performance:

Performance
===========

In order to speed up the simulations we here demonstrate how to run MARVELsim on a High Performace Computing (HPC) facility. We note that PyEchelle (and thus MARVELsim) is specifically designed to run with CUDA NVIDIA hardware which typically is available for GPU nodes on most computing clusters. Since PyEchelle is the natural bottleneck w.r.t. computation time, running large simulations using CUDA, or running the simulation on multiple CPUs, decrease the computational time dramatically. We specifically refer to `PyEchelle's documentation on performance <https://stuermer.gitlab.io/pyechelle/benchmark.html>`_ for more information. 

(and standard CPU parallelisation)


Job script for calibration mode
-------------------------------

We will use the Vlaams Supercomputing Centre (VSC) as example on how to create a job script that runs with GPUs. We likewise provide two job script examples in the ``MARVELsim/examples/HPC`` folder for your convenience. The first (``run_calibs.pbs``) shows a typical job script for running a full set of calibrated data:

.. code-block:: shell

   #!/bin/bash

   #PBS -N output
   #PBS -A <account_name>
   #PBS -l nodes=1:ppn=36:gpus=4:skylake
   #PBS -l partition=gpu
   #PBS -l pmem=2gb
   #PBS -l walltime=03:00:00

   cd $PBS_O_WORKDIR

   POETRY=/data/leuven/341/vsc34166/poetry/virtualenvs/marvelsim-3qcQCF7a-py3.8/bin/activate
   export POETRY
   SIMDIR=$VSC_DATA/MARVELsim/marvelsim
   export SIMDIR
   OUTDIR=/scratch/leuven/341/vsc34166/marvelsim/sim1/output
   export OUTDIR

   # Activate poetry shell
   source $POETRY
   
   # Run MARVELsim
   python $SIMDIR/marvelsim.py --calibs --cuda --zip -o $OUTDIR

Illustrated here we request a single node with 4 GPUs each using 9 CPU slaves (hence 36 in total) to execute the job. We request 2 GB of memory RAM to be on the safe side since a single 10,560 x 10,560 pixel full frame image occupy 851 Mb. Notice that it is possible to compress each image on the fly by enabling the flag ``zip`` as done in this example. Typical deflation rates per image are around 80%, hence, it is highly recommended to invoke this flag for faster data transfer after end job. For the job script show above the total run time (walltime) was 2 hours and 40 minutes.

Job script for science mode
---------------------------

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

   POETRY=/data/leuven/341/vsc34166/poetry/virtualenvs/marvelsim-3qcQCF7a-py3.8/bin/activate
   export POETRY
   SIMDIR=$VSC_DATA/MARVELsim/marvelsim
   export SIMDIR
   OUTDIR=/scratch/leuven/341/vsc34166/marvelsim/sim1/output
   export OUTDIR

   # Run MARVELsim shell
   source $POETRY
   python $SIMDIR/marvelsim.py --time 900 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --data rv_data.txt --cuda --zip -o $OUTPUT

Akin to the previous job script we here use the same computational resources, however, with the exception of increasing the walltime. Notice that adding more nodes will not speed up the computations, however, some cluster do provide more GPUs which will decrease the run time. We recommend to debug and test the computational resources needed for your jobs adding ``#PBS -l qos=debugging`` to the PSB details in the scripts shown above.  


Step-by-step guide
------------------

In order to run a simulation on the HPC the following needs to be secured:

 - Use the workflow script
 - Change the name of the data.txt file used to run with worker
 - Change the information within ``run_science_pyechelle.pbs`` and ``run_science_pyxel.pbs``
 - Change the ``ouput_folder`` within the file ``inputfiles/inputfile_marvelsim.yaml``
 - If saving the files to VSC scratch (as recommended) never use ``$VSC_SCRATCH`` but use e.g. ``/scratch/leuven/341/vsc34166/marvelsim/sim0``


VSC information
---------------

We note that the aboved resources w.r.t. skylake GPU nodes are the maximum and, hence, the computation times stated above using the VSC are at their minimum.

To get started using the VSC infrastrutrue we recommend reading:
  - `Genius quickstart guide <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/genius_quick_start.html#submit-to-genius-gpu-node>`_
  - `Genius hardware <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/tier2_hardware/genius_hardware.html>`_

