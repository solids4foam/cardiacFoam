#!/bin/bash

cd "$1"
source /Volumes/OpenFOAM-v2412/etc/bashrc

echo "🧹 Cleaning case"
./Allclean

echo "🚀 Running single Cell Solver"
./Allrun




