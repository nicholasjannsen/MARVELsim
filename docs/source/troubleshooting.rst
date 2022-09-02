Troubleshooting
===============

Please report bugs and issue via a GitHub issues through the `MARVELsim repository <https://github.com/nicholasjannsen/MARVELsim>`_.

Cluster complications
---------------------

- When saving files on a computing cluster you may experience that the symbolic link like ``$DATA``, ``$SCRATCH``, etc. are not globally accessable to the node-cores for which the simulations are taking place, even if these paths are being exported as part of your job scripts. E.g. on the VSC cluster the **scratch** file location (called ``VSC_SCRATCH``) are not recognized on the compute nodes and hence here we need to use the abosule path ``/scratch/path/to/marvelsim/output``. 
