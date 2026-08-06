import pyvista as pv
mesh = pv.read('exported_combined_purkinje/combined_step_050.vtm')
heart = mesh[0]
purkinje = mesh[1]

# Rename Vm_V to Vm so they align perfectly
if 'Vm_V' in purkinje.point_data:
    purkinje.point_data['Vm'] = purkinje.point_data.pop('Vm_V')

# Add region id
import numpy as np
heart.point_data["RegionId"] = np.zeros(heart.n_points, dtype=np.int32)
purkinje.point_data["RegionId"] = np.ones(purkinje.n_points, dtype=np.int32)

merged = heart.merge(purkinje)
print('Merged type:', type(merged))
print('Merged point data:', merged.point_data.keys())

merged.save('test_merged.vtp')
print('Saved test_merged.vtp')
