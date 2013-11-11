#!/bin/bash

L=$1

if [ "X$L" = "X" ]; then
	echo "Must specify L"
	exit
fi

IFILE="R/run_model2d.R"
OFILE="output/L$L.out"

# read in profile
. ~/.bash_profile

# move to root directory
cd ..

# execute R script
~/soft/lib64/R/bin/R -q --no-save --args "$L" < $IFILE > $OFILE 2>&1
