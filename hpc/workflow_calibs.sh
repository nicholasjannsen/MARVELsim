#!/bin/bash

# Clean and load modules
module purge
module restore plato
module load worker

# Summit jobs as a workflow
workflow1=$(qsub run_calibs_pyechelle.pbs)

wsub -W depend=afterok:$workflow1 -batch run_science_pyxel.pbs -data data.txt
