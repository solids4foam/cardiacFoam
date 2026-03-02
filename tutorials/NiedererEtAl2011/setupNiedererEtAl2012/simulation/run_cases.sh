#!/bin/bash

cd "$1"
source /Volumes/OpenFOAM-v2412/etc/bashrc

echo "🧹 Cleaning case"
./Allclean

echo "🚀 Running parallel"
./Allrun parallel

echo "✍️ Touching case.foam for paraview"
touch case.foam


