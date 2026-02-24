"""Default paths and parameters for scar tagging."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES_ORGANIZE = ROOT.parent / "filesOrganize"
OUTPUTS = FILES_ORGANIZE / "cardiac_preproc" / "outputs"

FULL_MESH = str(OUTPUTS / "purkinjeLayer_Diffusivity_IDsGlobal.vtk")
SELECTION = str(OUTPUTS / "scar_tissue_region.vtu")
OUTPUT = str(OUTPUTS / "purkinjeLayer_Diffusivity_scar.vtk")

# Diffusivity scaling for scar cells.
# mode: "constant" applies DIFFUSIVITY_SCAR_MULTIPLIER to scar cells.
DIFFUSIVITY_MODE = "constant"
DIFFUSIVITY_SCAR_MULTIPLIER = 0.0
