#!/bin/bash

# Parsed arguments
data=$1

# Clean and load modules
module purge
module restore marvelsim
module load worker

# Summit jobs as a workflow
workflow1=$(qsub run_science_pyechelle.pbs)
wsub -W depend=afterok:$workflow1 -batch run_science_pyxel.pbs -data $data
