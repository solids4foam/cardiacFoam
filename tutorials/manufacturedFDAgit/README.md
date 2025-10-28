The `monodomain_fv_cell_centred.cpp` program is a basic explicit, cell-centred
 implementation, which can be compiled and run with:

```
g++ -O2 -std=c++17 monodomain_fv_cell_centred.cpp -o monodomain_fv_cell_centred && ./monodomain_fv_cell_centred
```

The output should show second order accuracy as below:

```
   N      dx        dt       steps      errV_inf      pV      errU1_inf     pU1      errU2_inf     pU2
  50 0.020000   0.000358       558      0.000015    0.000      0.000080    0.000      0.000038    0.000
 100 0.010000   0.000090      2230      0.000004    2.000      0.000020    1.997      0.000010    1.999
 200 0.005000   0.000022      8917      0.000001    2.000      0.000005    1.999      0.000002    2.000
 400 0.002500   0.000006     35666      0.000000    2.000      0.000001    2.000      0.000001    2.000

Notes:
 - Cell-centred FV grid: x_i = (i+0.5)*dx on [0,1].
 - Explicit Euler in time with dt ~ dx^2 for diffusion stability.
```