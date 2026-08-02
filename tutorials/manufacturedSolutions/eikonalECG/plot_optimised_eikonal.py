import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv('setup/results/scheme_study_optimised.csv')

# Extract values
dx = df['dx'].values
at_l2 = df['activationTime_L2'].values
ecg_l2 = df['ecg_L2'].values

# Calculate theoretical slopes matching the first point
o1 = at_l2[0] * (dx / dx[0])**1.0
o2 = ecg_l2[0] * (dx / dx[0])**2.0

plt.figure(figsize=(8, 6))
plt.loglog(dx, at_l2, 'o-', linewidth=2, markersize=8, label='Activation Time $L_2$', color='blue')
plt.loglog(dx, ecg_l2, 's-', linewidth=2, markersize=8, label='Pseudo-ECG $L_2$', color='green')

# Reference lines
plt.loglog(dx, o1, 'k--', label='$\mathcal{O}(h^1)$ reference')
plt.loglog(dx, o2, 'k:', label='$\mathcal{O}(h^2)$ reference')

plt.xlabel('Mesh Spacing $\Delta x$ ($h$)')
plt.ylabel('Relative $L_2$ Error')
plt.title('Eikonal Convergence on Optimised Frontal Mesh')
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.tight_layout()

out_path = '/Users/simaocastro/.gemini/antigravity/brain/18d12a46-f88c-400c-8c11-37df34b853eb/eikonal_optimised_convergence.png'
plt.savefig(out_path, dpi=300)
print(f"Saved plot to {out_path}")
