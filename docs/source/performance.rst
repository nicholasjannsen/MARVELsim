Performance
===========

In order to speed up the simulations we here demonstrate how to run MARVELsim on a High Performace Computing (HPC) facility. We note that PyEchelle (and thus MARVELsim) is specifically designed to run with CUDA NVIDIA hardware which typically is available for GPU nodes on most computing clusters. Since PyEchelle is the natural bottleneck w.r.t. computation time, running large simulations using CUDA, or running the simulation on multiple CPUs, decrease the computational time dramatically. We specifically refer to `PyEchelle's documentation on performance <https://stuermer.gitlab.io/pyechelle/benchmark.html>`_ for more information. 

We will use the Vlaams Supercomputing Centre (VSC) as example on how to create a job script that runs with GPUs. We likewise provide two job script examples in the ``MARVELsim/examples/HPC`` folder for your convenience. For each job script we specify the resources (see ``#PBS ..``), next we export all paths needed to be recognized globally from any compute node, and lastly we activate our python environment and launch MARVELsim.  

.. _performance_calib_mode:

Job script - Calibration mode
-----------------------------

The first example (``run_calibs.pbs``) shows a typical job script for running a full set of calibrated data:

.. code-block:: shell

   #!/bin/bash

   #PBS -N output
   #PBS -A <account_name>
   #PBS -l nodes=1:ppn=36:gpus=4:skylake
   #PBS -l partition=gpu
   #PBS -l pmem=2gb
   #PBS -l walltime=03:00:00

   cd $PBS_O_WORKDIR

   POETRY=$VSC_DATA/poetry/virtualenvs/marvelsim-3qcQCF7a-py3.8/bin/activate
   export POETRY
   SIMDIR=$VSC_DATA/MARVELsim/marvelsim
   export SIMDIR
   OUTDIR=/scratch/leuven/341/vsc34166/marvelsim/output
   export OUTDIR

   # Activate poetry shell
   source $POETRY
   
   # Run calibration mode
   python $SIMDIR/marvelsim.py --calibs --cuda --zip -o $OUTDIR

Illustrated here activate the flag ``--cuda`` and request a single node with 4 GPUs each using 9 CPU slaves (hence 36 in total) to execute the job. We request 2 GB of memory RAM to be on the safe side since a single 10,560 x 10,560 pixel full frame image occupy 851 Mb. In order to activate your Poetry shell the absolute path needs to be exported globally. If ``poetry shell`` is activated you can simply type ``which python`` to get the absolute path needed to add to the above variable ``POETRY``. We here save the output data to the **scratch** file location in order avoid overflowing our memory storage on the **data** storage.

Similar to the example given in the :ref:`tutorial <tutorial_calibration>` we use the off-the-shelf method of producing a set of calibrated data simply by invoking the flag ``--calibs``. Notice that it is possible to compress each image on the fly by enabling the flag ``zip`` as done in this example. Typical deflation rates per image are around 80%, hence, it is highly recommended to invoke this flag for faster data transfer after end job. For the job script shown above the total run time (a.k.a. walltime) was 2 hours and 40 minutes.


.. _performance_science_mode:

Job script - Science mode
-------------------------

The following example (``run_science.pbs``) shows a job script for running 300 stellar spectra using a generated RV time series called ``rv_data.txt``:

.. code-block:: shell

   #!/bin/bash

   #PBS -N output
   #PBS -A <account/project>
   #PBS -l nodes=1:ppn=36:gpus=4:skylake
   #PBS -l partition=gpu
   #PBS -l pmem=2gb
   #PBS -l walltime=10:00:00

   cd $PBS_O_WORKDIR

   POETRY=$VSC_DATA/poetry/virtualenvs/marvelsim-3qcQCF7a-py3.8/bin/activate
   export POETRY
   SIMDIR=$VSC_DATA/MARVELsim/marvelsim
   export SIMDIR
   OUTDIR=/scratch/leuven/341/vsc34166/marvelsim/output
   export OUTDIR

   # Activate poetry shell
   source $POETRY
   
   # Run science mode
   python $SIMDIR/marvelsim.py --science --time 900 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --data rv_data.txt --cuda --zip -o $OUTPUT

Akin to the previous job script we here use the same computational resources, however, with the exception of increasing the walltime. Notice that adding more nodes will not speed up the computations, however, some cluster do provide more GPUs which will decrease the run time. We recommend to debug and test the computational resources needed for your jobs adding ``#PBS -l qos=debugging`` to the PSB details in the scripts shown above and run a single simulation.  

Workflow - Science mode
-----------------------

While PyEchelle runs very efficiently using GPUs, the code parallisation for Pyxel is most efficient while using normal CPUs. Thus to get the minimum walltime for your simulations we recomment to use a workflow, i.e. a step-wise execution of two or more codes. Hence, in the following we show a workflow script (``worflow_science.sh``) that first executes the PyEchelle simulations using GPUs, and only when finished successfully, the job script will continue to launch the Pyxel simulations using CPUs:

.. code-block:: shell

   #!/bin/bash
		
   # First summit PyEchlle job
   workflow=$(qsub run_science_pyechelle.pbs)

   # When finished successfully summit Pyxel job
   wsub -W depend=afterok:$workflow -master -batch run_science_pyxel.pbs -data rv_data.txt

Like before we here used the standard Torque schedular command ``qsub`` to summit the PyEchelle job. The Pyxel job is submitted using the popular ``worker`` framework. By default worker use one node-core to schedule the simulation, however, as we only have a smaller amount of jobs (300 in total) we can overwrite this behavior and tell worker to use all node-cores for the computation. This is simply done by using the flag ``-master``. Worker will automatically parameterise the ``rv_data.txt`` file for which we use the index and the RV amplitude from (see the output of the :ref:`RV generator <tutorial_rv_script`).


  

Step-by-step guide
------------------

In order to run a simulation on any cluster the following needs to be secured:

 - Copy the job script and RV data to your cluster
 - If generating science simulations consider using workflows (see job scripts below)
 - Adjust the job script details: resources, paths/names, input parameters, etc.
 - Adjust output destination of your simulation. Notice that ``OUTDIR`` within your job script and the ``ouput_folder`` within the input file ``inputfiles/inputfile_marvelsim.yaml`` need to match and can only be absolute paths (hence do not use symbolic links like ``$DATA``) 


VSC information
---------------

- We note that the aboved resources w.r.t. skylake GPU nodes are the maximum and, hence, the computation times stated above using the VSC are at their minimum.
- To get started using the VSC infrastrutrue we recommend reading:
  - `Genius quickstart guide <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/genius_quick_start.html#submit-to-genius-gpu-node>`_
  - `Genius hardware <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/tier2_hardware/genius_hardware.html>`_

