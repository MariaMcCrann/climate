#!/bin/bash

while read L; do
	R -q --no-save --args "$L" <R/run_model2d.R
done <Ls
