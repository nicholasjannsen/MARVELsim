README
======

This note serves a some advices and important information in order to run simulations on the HPC successfully.

In order to run a simulation on the HPC the following needs to be secured:

 - Use the workflows script
 - Change the name of the data.txt file used to run with worker
 - Change the information within ``run_science_pyechelle.pbs`` and ``run_science_pyxel.pbs``
 - Change the ``ouput_folder`` within the file ``inputfiles/inputfile_marvelsim.yaml``
 - If saving the files to VSC scratch (as recommended) never use ``$VSC_SCRATCH`` but use e.g. ``/scratch/leuven/341/vsc34166/marvelsim/sim0``
   
Running with the default workflows given in the repo the computational resources are the following:

Science:
 - Time: 
