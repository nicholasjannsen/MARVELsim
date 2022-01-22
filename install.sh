#/usr/bin/env bash

export MARVELSIM=$PWD
source $MARVELSIM/marvelsim/bin/activate

pip install numpy==1.21
pip install numba==0.55
pip install colorama
pip install pyechelle
pip install pyxel-sim[all]==11.5
