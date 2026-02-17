"""Default paths for VTK FIELD conversion step."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = ROOT / "outputs"

INPUT = str(OUTPUTS / "purkinjeLayer_Diffusivity_scar.vtk")
OUTPUT = str(OUTPUTS / "purkinjeLayer_Diffusivity_scar_VtkUnstructured.vtk")
