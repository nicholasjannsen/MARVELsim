Troubleshooting
===============

Please report bugs and issue via a GitHub issues through the `MARVELsim repository <https://github.com/nicholasjannsen/MARVELsim>`_.

Pyxel's output path
-------------------

Note that while working in a python prompt/environment, Pyxel currently has a bug for overwriting the output file directory stated in the input YAML file. This means that the user needs to manually specify the output path within ``MARVELsim/inputfiles/inputfile_marvel.yaml`` (see ``exposure:outputs:output_folder``). Otherwise the default ``MARVELsim/output`` folder will be used. If you try to save data 

Cluster complications
---------------------

- When saving files on a computing cluster you may experience that the symbolic link like ``$DATA``, ``$SCRATCH``, etc. are not globally accessable to the node-cores for which the simulations are taking place, even if these paths are being exported as part of your job scripts. E.g. on the VSC cluster the **scratch** file location (called ``VSC_SCRATCH``) are not recognized on the compute nodes and hence here we need to use the abosule path ``/scratch/path/to/marvelsim/output``. 
