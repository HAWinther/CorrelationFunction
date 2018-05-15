#!/bin/bash

#=======================================================
# Testrun of code. Expected output is given in [corr.out]
# Assuming OpenMP version with weights
#=======================================================

codedir="../"
thisdir=$(pwd)

# Compile
cd $codedir
make clean; make
cp CUTER $thisdir/
cd $thisdir

# Run
galaxyfile="gal.dat"
randomfile="ran.dat"
outfile="corr.out"
nbins="32"
rmax="192.0"
OmegaM="0.3"

./CUTER $galaxyfile $randomfile $outfile $nbins $rmax $OmegaM

