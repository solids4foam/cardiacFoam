# Extracting 3D Time-Series Data (Probes)

When running 3D myocardium simulations, `cardiacFoam` automatically outputs the requested ionic variables (like `Vm`, `Cai`, `CaSR`) as native OpenFOAM `volScalarField`s. 

To prevent generating massive amounts of disk data while still getting high-resolution time-series curves (e.g., to plot a smooth Action Potential or Calcium transient), you should use the native OpenFOAM `probes` function. 

## How to use `probes`
You do not need to write any C++ code. You simply add a `probes` block to your `system/controlDict` (or include it as a separate file, like is done in `tutorials/NiedererEtAl2011/tissueNiedererEtAl2011/system/Niedererpoints`).

Here is a complete example of how to extract data for two specific cells every single timestep:

```cpp
functions
{
    myProbes
    {
        type            probes;
        libs            ("libsampling.so");
        
        // Write frequency (timeStep 1 = every single time step)
        writeControl    timeStep;
        writeInterval   1;
        
        // The specific fields you want to track over time
        fields
        (
            Vm
            Cai
        );
        
        // The (x y z) coordinates of the cells you want to extract
        probeLocations
        (
            (0.01  0.0  0.0)    // e.g. Endocardial cell
            (0.05  0.0  0.0)    // e.g. Mid-myocardial cell
        );
    }
}
```

When you run your solver, OpenFOAM will create a `postProcessing/myProbes/0/` directory containing beautifully formatted `.dat` files for each field, which you can immediately drop into Python, Excel, or gnuplot.
