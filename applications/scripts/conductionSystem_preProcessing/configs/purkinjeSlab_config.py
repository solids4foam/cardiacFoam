"""
config.py

Configuration for Purkinje slab tagging.
"""

# --- Purkinje layer defaults (layer mode) ---

    #Layer thickness as fraction of transmural depth( 10% transmural depth starting in subendocardium)
PURKINJE_TRANSMURAL_MIN = 0.0
PURKINJE_TRANSMURAL_MAX = 0.1

    #Tags for LV and RV in uvc_intraventricular field (Check fields list to confirm)
LV_INTRAVENTRICULAR_TAG = -1
RV_INTRAVENTRICULAR_TAG =  1 



    #Named for the cell data field tagging Purkinje elements for the slab approach
PURKINJE_FIELD_NAME = "purkinjeLayer"

# --- Diffusivity scaling for Purkinje slab ---
PURKINJE_DIFFUSION_MULTIPLIER = 4.0

# --- Optional wall/boundary filtering ---
# Use a tag from cell_data (e.g. "endo_surface_lv") or field_data + gmsh:physical.
# Keep as None to disable wall restriction.
WALL_TAG_NAME = None
WALL_TAG_VALUE = None

# --- Optional inflate-from-point geometry ---
# Seed point "x,y,z" to grow region until shared boundary (same rule as tag_inside_shared_boundary).
# Keep as None to disable.
INFLATE_SEED_POINT = None
INFLATE_INSIDE_TAG_NAME = "inside_shared_boundary"

# --- Default input/output paths ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "inputs" / "meshes"
OUTPUTS = ROOT / "outputs"

INPUT = str(INPUTS / "ASCIIlegacy.vtk")
OUTPUT = str(OUTPUTS / "purkinjeLayer.vtk")
