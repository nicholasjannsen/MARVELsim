.. _performance:

Performance
===========

In order to speed up the simulations we here demonstrate how to run MARVELsim on a High Performace Computing (HPC) facility. Generally, since PyEchelle is the natural bottleneck w.r.t. computation time, the very user friendly PyEchelle interface to run on GPUs or normal CPUs is beneficially used in MARVELsim. As this dramatically decrease the run time, we strongly recommend to checkout `PyEchelle's documentation on performance <https://stuermer.gitlab.io/pyechelle/benchmark.html>`_ for a deeper understanding of what's going on under the hood in what follows. We will use the Vlaams Supercomputing Centre (VSC) as example on how run the calibration- and science mode of MARVELsim.

Workflows
---------

If available PyEchelle is extremely efficient to run with CUDA on NVIDIA hardware which typically is available for GPU nodes on most computing clusters. On the other hand Pyxel is not developed for the usage of GPUs but rather for normal CPU prallelisation. Thus, to not waste unesseary computional resources, we will in the following show how to run a so-called *workflow*; that is, to summit a combined script that first runs software 1 (i.e. PyEchelle on GPUs), and only when this finish succesfully, then run software 2 (i.e. Pyxel on CPUs) that has a input dependency from software 1 (i.e. the CCD full-frame images). We conveniece we provide a ready-to-go script to be executed on the VSC:

.. code-block:: shell

   #!/bin/bash                                                                                                                                 

   # Clean and load modules                                                                                                                    
   module purge
   module restore plato
   module load worker

   # Summit jobs as a workflow                                                                                                                 
   workflow1=$(qsub run_science_pyechelle.pbs)
   wsub -W depend=afterok:$workflow1 -batch run_science_pyxel.pbs -data data_200kms.txt

Currently, we only provide a workflow script (``MARVELsim/hpc/workflow_science.sh``) for the science mode. The important details here are the two job scripts called ``run_science_pyechelle.pbs`` and ``run_science_pyxel.pbs`` which each will invoke MARVELsim to run each software individually. We explain the details of these in the follwoing. 


Job script for science mode
---------------------------

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
   python simulator-marvel.py --time 900 --mag 10.0 --teff 5800 --logg 4.5 --z 0.0 --alpha 0.0 --data rv_data.txt --cuda -o $SIMDIR/output

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
		
Job script for calibration mode
-------------------------------

We likewise provide two job script to run a set of calibrated images called ``run_calibs.pbs`` on the VSC:

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


Compared to the science mode we haven't made an effort to split up the computation between previous job scripts we here use the same computational resources


VSC information
---------------

Notice that adding more nodes will not speed up the computations, however, some cluster do provide more GPUs which will decrease the run time.
We recommend to debug and test the computational resources needed for your jobs adding ``#PBS -l qos=debugging`` to the PSB details in the scripts shown above.  

We note that the aboved resources w.r.t. skylake GPU nodes are the maximum and, hence, the computation times stated above using the VSC are at their minimum.

To get started using the VSC infrastrutrue we recommend reading:
  - `Genius quickstart guide <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/genius_quick_start.html#submit-to-genius-gpu-node>`_
  - `Genius hardware <https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/leuven/tier2_hardware/genius_hardware.html>`_

