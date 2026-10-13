#!/usr/bin/env python3
"""Independent 1D reciprocal-current/area and Wilson-reference check."""
from pathlib import Path
import re

def vectors(name):
    text=Path('0',name).read_text().split('internalField',1)[1].split('boundaryField',1)[0]
    return [tuple(map(float,x.split())) for x in re.findall(r'\(([^()]+)\)',text)]

b=vectors('leadField_bipolar'); w=vectors('leadField_V1')
assert len(b)==len(w)==60
# In the bulk of the 1m heart column: unit reciprocal current / area .01,
# total Gi+Ge=4 S/m, and Gi=2 S/m. Lx=Gi/(area*(Gi+Ge))=50 m^-2.
for i in range(20,40):
    assert abs(b[i][0]-50)<1e-9 and max(abs(v) for v in b[i][1:])<1e-12, (i,b[i])
# All three Wilson reference cells lie to the right of the entire heart;
# their summed weight is -1, so the myocardial sensitivity matches bipolar.
for i in range(20,40):
    assert max(abs(a-c) for a,c in zip(b[i],w[i]))<1e-10*max(1,max(map(abs,b[i]))), (i,b[i],w[i])
print('PASS: reciprocal-current normalization, polarity and Wilson weights')
