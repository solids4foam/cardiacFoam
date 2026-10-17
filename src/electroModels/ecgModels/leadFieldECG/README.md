# leadFieldECG

The ECG of a heart+torso passive conductor from precomputed lead fields
(`buildLeadFields`): each channel is the sum, over the heart cells, of the
electrode's lead field times the heart source of the current `Vm`.

Electrodes lie in the torso; intracardiac electrodes (in a myocardium cell)
are not supported and stop the run with a fatal error.
