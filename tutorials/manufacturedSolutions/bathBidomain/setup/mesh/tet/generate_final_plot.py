import os
import vtk
import numpy as np
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy

def get_mesh_metrics(path, max_cells=1000000):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(path)
    reader.Update()
    grid = reader.GetOutput()
    
    points = vtk_to_numpy(grid.GetPoints().GetData())
    cell_types = vtk_to_numpy(grid.GetCellTypes())
    cell_array = grid.GetCells()
    conn = vtk_to_numpy(cell_array.GetConnectivityArray())
    offsets = vtk_to_numpy(cell_array.GetOffsetsArray())
    
    sizes = np.diff(offsets)
    is_tet = (cell_types == 10) & (sizes == 4)
    starts = offsets[:-1][is_tet]
    idx = starts[:, None] + np.arange(4)[None, :]
    tets = conn[idx][:max_cells]
    
    cc = points[tets].mean(axis=1)
    combos = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])
    faces = tets[:, combos].reshape(-1, 3)
    owner = np.repeat(np.arange(len(tets)), 4)
    
    key = np.sort(faces, axis=1)
    order = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    key_s, owner_s, faces_s = key[order], owner[order], faces[order]
    
    same = np.all(key_s[:-1] == key_s[1:], axis=1)
    first = np.flatnonzero(same)
    own = owner_s[first]
    nei = owner_s[first + 1]
    fverts = faces_s[first]
    
    p0, p1, p2 = (points[fverts[:, i]] for i in range(3))
    fctr = (p0 + p1 + p2) / 3.0
    sf = 0.5 * np.cross(p1 - p0, p2 - p0)
    magsf = np.linalg.norm(sf, axis=1)
    d = cc[nei] - cc[own]
    magd = np.linalg.norm(d, axis=1)
    
    ok = (magsf > 0) & (magd > 0)
    sf, magsf, d, magd = sf[ok], magsf[ok], d[ok], magd[ok]
    fctr, own_c = fctr[ok], cc[own][ok]
    
    cosang = np.abs(np.einsum("ij,ij->i", d, sf)) / (magd * magsf)
    nonorth = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
    
    dn = np.einsum("ij,ij->i", d, sf)
    safe = np.abs(dn) > 1e-300
    cpf = fctr - own_c
    t = np.zeros(len(d))
    t[safe] = np.einsum("ij,ij->i", cpf[safe], sf[safe]) / dn[safe]
    pierce = own_c + t[:, None] * d
    skew = np.linalg.norm(fctr - pierce, axis=1) / magd
        
    return nonorth, skew

print("Reading Original...")
no_orig, sk_orig = get_mesh_metrics("original_N20.vtk")
print("Reading Optimized...")
no_opt, sk_opt = get_mesh_metrics("optimized_N20.vtk")

print("Reading Pig1...")
no_pig, sk_pig = get_mesh_metrics("/Users/simaocastro/pigData/pig1/mesh.vtk", max_cells=1000000)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

labels = ['Pig Anatomy\n(Target)', 'Original\nTet MMS', 'Optimized\nTet MMS']
colors = ['#1f77b4', '#d62728', '#2ca02c'] # Blue, Red, Green

def plot_violins(ax, data, title, ylabel, ymax, fmt='.1f'):
    parts = ax.violinplot(data, showmeans=False, showextrema=False)
    for pc, color in zip(parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')
        
    for i, d in enumerate(data):
        mean_val = np.mean(d)
        ax.plot([i+0.8, i+1.2], [mean_val, mean_val], color='black', lw=2)
        
        # Shift text slightly up by adding a tiny fraction of the axis range so it floats above the line
        y_offset = ymax * 0.015
        ax.text(i+1.0, mean_val + y_offset, f'{mean_val:{fmt}}', va='bottom', ha='center', fontsize=11, fontweight='bold')
        
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(labels, fontsize=12)
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.set_ylim(0, ymax)

plot_violins(ax1, [no_pig, no_orig, no_opt], 'Non-Orthogonality Distribution', 'Angle (degrees)', 90, fmt='.1f')

# Create custom legend
import matplotlib.patches as mpatches
legend_patches = [mpatches.Patch(color=c, label=l.replace('\n', ' ')) for c, l in zip(colors, labels)]
fig.legend(handles=legend_patches, loc='lower center', ncol=3, fontsize=12, bbox_to_anchor=(0.5, -0.05))

plot_violins(ax2, [sk_pig, sk_orig, sk_opt], 'Skewness Distribution', 'Skewness', 2.5, fmt='.3f')

plt.tight_layout()
plt.savefig("final_mesh_quality_plot.png", dpi=300, bbox_inches='tight')
print("Plot saved to final_mesh_quality_plot.png")
