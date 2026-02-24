"""Default paths and parameters for scar tagging."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES_ORGANIZE = ROOT / "filesOrganize"
CARDIAC_OUTPUTS = FILES_ORGANIZE / "cardiac_preproc" / "outputs"

FULL_MESH = str(CARDIAC_OUTPUTS / "purkinjeLayer_Diffusivity_IDsGlobal.vtk")
SELECTION = str(CARDIAC_OUTPUTS / "scar_tissue_region.vtu")
OUTPUT = str(CARDIAC_OUTPUTS / "purkinjeLayer_Diffusivity_scar.vtk")

DIFFUSIVITY_MODE = "constant"
DIFFUSIVITY_SCAR_MULTIPLIER = 0.0
