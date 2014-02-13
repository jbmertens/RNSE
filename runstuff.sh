#!/bin/bash

WDIR=data/$1

make clean
make fast=1
rm -r $WDIR

./rnse -a 0.90 -t 64 -c 0.1 -o $WDIR

export DISPLAY=:1
cp gifgen.math $WDIR/gifgen.math
cd $WDIR
math < gifgen.math

cd images
avconv -r 10 -f image2 -i slices-%d.png -b 4096k viscev.mp4

echo "done"
