#!/usr/bin/env bash

# Account details
VSC="$PWD/server"
X=341
Y=66

# Provide directory structure
if [ ! -d "$VSC" ]; then
	echo "Creating VSC file structure locally"
	mkdir $VSC
	mkdir $VSC/home
	mkdir $VSC/data
	mkdir $VSC/scratch
fi

# If no arguments are given write usage message
if [ -z "$1" ]; then
    echo "Usage: server-sshfs <option>"
	echo "       mount    :   <option> = m"
	echo "       dismount :   <option> = d"
	exit 1

else

    # Mount VSC account
    if [ $1 = "m" ]; then
        echo "Mounting VSC account"
		sshfs vsc$X$Y@login.hpc.kuleuven.be:/user/leuven/$X/vsc$X$Y $VSC/home
		sshfs vsc$X$Y@login.hpc.kuleuven.be:/data/leuven/$X/vsc$X$Y $VSC/data
		sshfs vsc$X$Y@login.hpc.kuleuven.be:/scratch/leuven/$X/vsc$X$Y $VSC/scratch
    fi

    # Dismount VSC account
    if [ $1 = "d" ]; then
        echo "Dismounting VSC account"
		fusermount -u $VSC/home
		fusermount -u $VSC/data
		fusermount -u $VSC/scratch
    fi
fi
