#!/bin/bash

# First summit PyEchlle job
workflow=$(qsub run_science_pyechelle.pbs)

# When finished successfully summit Pyxel job
wsub -W depend=afterok:$workflow -master -batch run_science_pyxel.pbs -data rv_data.txt
