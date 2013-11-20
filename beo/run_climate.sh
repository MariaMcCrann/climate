#!/bin/bash

L=$1
WC=$2

if [ "X$L" = "X" ]; then
	echo "Must specify L"
	exit
fi

if [ "X$WC" = "X" ]; then
	echo "Warning: Using default of summer temperature data."
	WC="ST"
fi

IFILE="R/run_model2d.R"
OFILE="output/L${L}_${WC}.out"

# read in profile
. ~/.bash_profile

# move to root directory
cd ..

# execute R script
~/soft/lib64/R/bin/R -q --no-save --args "$L $WC" < $IFILE > $OFILE 2>&1
