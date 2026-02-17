"""
config.py

Configuration for diffusivity tensors.
"""

# --- Diffusivities (ventricular tissue) in SI units: S/m ---
Ventricular_DF = 0.1143  # Along fiber direction
Ventricular_DS = 0.052  # Along sheet direction
Ventricular_DN = 0.052  # Normal to sheet direction

# --- Default input/output paths ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"

INPUT = str("ASCIIlegacy.vtk")
OUTPUT = str(OUTPUTS / "ASCIIlegacy_DIffusionTensor.vtk")
z