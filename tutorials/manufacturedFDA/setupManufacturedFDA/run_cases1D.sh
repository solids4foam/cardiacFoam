#!/bin/bash

cd "$1"
source /Volumes/OpenFOAM-v2412/etc/bashrc

echo "🧹 Cleaning case"
./Allclean

echo "🚀 Running "
./Allrun parallel

echo "✍️ extracting errors in main"



