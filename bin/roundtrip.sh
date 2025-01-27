#!/bin/bash

set -x
./jcodec_main 0 0 $1 $2 $3 $4
oname=$(basename $1 | sed "s/\.ppm//")
ppmname=${oname}_$4_deco.ppm
oname=${oname}_$4.vtp
./jcodec_main 1 0 $oname ppm
./ssim.py ${ppmname} $1
