"""Default paths and material parameters for diffusivity tagging."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES_ORGANIZE = ROOT / "filesOrganize"
CARDIAC_INPUTS = FILES_ORGANIZE / "cardiac_preproc" / "inputs" / "meshes"
CARDIAC_OUTPUTS = FILES_ORGANIZE / "cardiac_preproc" / "outputs"

# Diffusivities (ventricular tissue) in SI units: S/m.
Ventricular_DF = 0.1143
Ventricular_DS = 0.052
Ventricular_DN = 0.016

INPUT = str(CARDIAC_INPUTS / "ASCIIlegacy.vtk")
OUTPUT = str(CARDIAC_OUTPUTS / "ASCIIlegacy_DIffusionTensor.vtk")
