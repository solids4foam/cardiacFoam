

## Post-processing and verification

Results can be post-processed to extract activation times and compared against
published benchmark data or other simulation codes.
The code presented in the setup runs several simulations and run the post-process
for mesh and time dependant activation times. 

### Plotting active tension outputs

The script `plotTensionVoltage` plots columns from one or two files using
header names (first row). It supports single-variable plots, scaling, and
overlaying a second file.

Examples:

```bash
./plotTensionVoltage activeTension.txt Ta
./plotTensionVoltage activeTension.txt Ta 1/60
./plotTensionVoltage activeTension.txt AV_Ca_i 1 AV_Tension 1
./plotTensionVoltage file1.txt AV_Ca_i 1 file2.txt Ca_i 1
```

To list available header names and their column indices:

```bash
./plotTensionVoltage --list activeTension.txt
```

---
