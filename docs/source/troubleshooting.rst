Troubleshooting
===============

Please report bugs and issue via a GitHub issues through the `MARVELsim repository <https://github.com/nicholasjannsen/MARVELsim>`_.

Pyxel's output path
-------------------

Note that while working in a python prompt/environment, Pyxel currently has a bug for overwriting the output file directory stated in the input YAML file. This means that the user needs to manually specify the output path within ``MARVELsim/inputfiles/inputfile_marvel.yaml`` (see ``exposure:outputs:output_folder``). Otherwise the default ``MARVELsim/output`` folder will be used. If you try to save data 
