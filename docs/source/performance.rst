Performance
===========

In order to speed up the simulations we here demonstrate how to run MARVELsim on a High Performace Computing (HPC) facility. We note that PyEchelle (and thus MARVELsim) is specifically designed to run with CUDA NVIDIA hardware which typically is available for GPU nodes on most computing clusters. Since PyEchelle is the natural bottleneck w.r.t. computation time, running large simulations using CUDA, or running the simulation on multiple CPUs, decrease the computational time dramatically. We specifically refer to `PyEchelle's documentation on performance <https://stuermer.gitlab.io/pyechelle/benchmark.html>`_ for more information. 

We will use the Vlaams Supercomputing Centre (VSC) as example on how to create a job script that runs with GPUs. We likewise provide two job script examples in the ``MARVELsim/examples/HPC`` folder for your convenience. For each job script we specify the resources (see ``#PBS ..``), next we export all paths needed to be recognized globally from any compute node, and lastly we activate our python environment and launch MARVELsim.  

.. admonition:: Step-by-step guide

   In order to run a simulation on any cluster the following needs to be secured:

   - If using science mode, generate a RV time series using ``rv-generator.py``
   ..
      - If using science mode, consider using a workflow (see job scripts below)
   - Copy one of the job script examples within ``MARVELsim/examples/HPC``
   - Adjust the job script details: resources, paths/names, input parameters, etc.
   - Adjust output destination of your simulation. Notice that ``OUTDIR`` within your job script and the ``ouput_folder`` within the input file ``inputfiles/inputfile_marvelsim.yaml`` need to match and can only be absolute paths (hence do not use symbolic links like ``$DATA``) 

     
.. _performance_calibs:

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


.. _performance_science:

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

Akin to the previous job script we here use the same computational resources, however, with the exception of increasing the walltime and the flag ``--science``. Notice that adding more nodes will not speed up the computations, however, some cluster do provide more GPUs which will decrease the run time. We recommend to debug and test the computational resources needed for your jobs adding ``#PBS -l qos=debugging`` to the PSB details in the scripts shown above and run a single simulation.  


..
   .. _performance_workflow:

   Workflow - Science mode
   -----------------------

   If available PyEchelle is extremely efficient to run with CUDA on NVIDIA hardware which typically is available for GPU nodes on most computing clusters. On the other hand Pyxel is not developed for the usage of GPUs but rather for normal CPU prallelisation. Thus, to not waste unesseary computional resources, we will in the following show how to run a so-called *workflow*; that is, to summit a combined script that first runs software 1 (i.e. PyEchelle on GPUs), and only when this finish succesfully, then run software 2 (i.e. Pyxel on CPUs) that has a input dependency from software 1 (i.e. the CCD full-frame spectra). For your conveniece we provide a ready-to-go script (``worflow_science.sh``) to be executed on the VSC:

   .. code-block:: shell

      #!/bin/bash

      # First summit PyEchlle job
      workflow=$(qsub run_science_pyechelle.pbs)

      # When finished successfully summit Pyxel job
      wsub -W depend=afterok:$workflow -master -batch run_science_pyxel.pbs -data rv_data.txt

   Like before we here used the standard Torque schedular command ``qsub`` to summit the PyEchelle job. The Pyxel job is submitted using the popular ``worker`` framework. By default worker use one node-core to schedule the simulation, however, as we only have a smaller amount of jobs (300 in total) we can overwrite this behavior and tell worker to use all node-cores for the computation. This is simply done by using the flag ``-master``. Worker will automatically parameterise the ``rv_data.txt`` file for which we use the index and the RV amplitude from (see the output of the :ref:`RV generator <tutorial_rv_script`).

   Currently, we only provide a workflow script (``MARVELsim/hpc/workflow_science.sh``) for the science mode. The important details here are the two job scripts called ``run_science_pyechelle.pbs`` and ``run_science_pyxel.pbs`` which each will invoke MARVELsim to run each software individually. We explain the details of these in the follwoing. 

   The following example (``run_science_pyechelle.pbs``) shows a job script for running 300 stellar spectra using a generated RV time series called ``rv_data.txt``:

   .. code-block:: shell

      #!/bin/bash

      #PBS -N output
      #PBS -A <account/project>
      #PBS -l nodes=1:ppn=36:gpus=4:skylake
      #PBS -l partition=gpu
      #PBS -l pmem=2gb
      #PBS -l walltime=10:00:00

      cd $PBS_O_WORKDIR

      PYTHONPATH=$VSC_DATA/MARVELsim/marvelsim/bin/python
      export PYTHONPATH
      SIMDIR=$VSC_DATA/MARVELsim
      export SIMDIR

      # Activate environment 
      source marvelsim/bin/activate

      # Run MARVELsim for PyEchelle only
      cd $SIMDIR
      python simulator-marvel.py --time 300 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --data rv_data.txt --cuda -o $SIMDIR/output

   Illustrated here we request a single node with 4 GPUs using each using 9 CPU claves (hence 36 in total) to execute the job. We request 2 GB of memory RAM to be on the safe side since a single 10,560 x 10,560 pixel full frame image occupy 851 Mb. The the run time (a.k.a. walltime) has here been timed to be around 10 hours.

   Next we call MARVELsim to invoke Pyxel only using the job script (``run_science_pyxel.pbs``):

   .. code-block:: shell

      #!/bin/bash                                                                                                                                 

      #PBS -N output                                                                                                                              
      #PBS -A <account>
      #PBS -l nodes=1:ppn=6:skylake                                                                                                               
      #PBS -l pmem=30gb                                                                                                                           
      #PBS -l walltime=04:00:00                                                                                                                   

      cd $PBS_O_WORKDIR

      PYTHONPATH=$VSC_DATA/MARVELsim/marvelsim/bin/python
      export PYTHONPATH
      SIMDIR=$VSC_DATA/MARVELsim
      export SIMDIR

      # Make sure to activate environment                                                                                                         
      source marvelsim/bin/activate

      # Run star spectrum                                                                                                                         
      cd $SIMDIR
      python simulator-marvel.py --time 900 --dex $index --zip -o $SIMDIR/output

   Seen here we only use 6 CPUs since Pyxel needs a very large amount of RAM memory for each image (of the order of 25 Gb), hence, using only 1 node we are limited here to 6 CPUs in order not to overflow the node memory. Notice that it is possible to compress each image on the fly by enabling the flag ``zip`` as done in this example. Typical deflation rates per image are around 80%, hence, it is highly recommended to invoke this flag for faster data transfer after end job. For the job script show above the total run time (walltime) was 3 hours. We further remark that Pyxel only needs the exposure time to apply CCD effects correctly which explains the absence of the stellar parameters. As shown from the workflow script above we used the popular *worker* framework to parallelise our simulations. Worker can immediately recognize the indices given in the first column of the RV data file ``rv_data.txt`` and used the ``$index`` parametrisation to automatically deligate the work to multiple CPU slaves.  


VSC information
---------------

- We note that the aboved resources w.r.t. skylake GPU nodes are the maximum and, hence, the computation times stated above using the VSC are at their minimum.
- To get started using the VSC infrastrutrue we recommend reading:
  - `Genius quickstart guide <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/genius_quick_start.html#submit-to-genius-gpu-node>`_
  - `Genius hardware <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/tier2_hardware/genius_hardware.html>`_

