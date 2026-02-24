"""Default paths for VTK FIELD conversion step."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FILES_ORGANIZE = ROOT.parent / "filesOrganize"
OUTPUTS = FILES_ORGANIZE / "cardiac_preproc" / "outputs"

INPUT = str(OUTPUTS / "purkinjeLayer_Diffusivity_scar.vtk")
OUTPUT = str(OUTPUTS / "purkinjeLayer_Diffusivity_scar_VtkUnstructured.vtk")
