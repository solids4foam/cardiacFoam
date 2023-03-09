# cardiacFoam
This repository contains *work in progress* on the implementation of cardiac electro-mechanical models in OpenFOAM.

# What does it contain?
Currently, the monodomain reaction-diffusion equation is implemented using the minimal Bueno-Orovio ionic model. This solver can be run on the included Niederer et al. (2011) benchmark case.

# How to install it?
The code requires OpenFOAM-v2012, although it should be relatively easy to port to other ESI, Foundation and foam-extend versions.
Once OpenFOAM-v2012 is sourced, you can compile the current toolbox using the included script:

    > ./Allwmake

You can then try the tutorial:

    > cd tutorials/NiedererBenchmarkSlab
    > ./Allrun

# Who to contact about it?
philip.cardiff@ucd.ie