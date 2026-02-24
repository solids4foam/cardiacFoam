# Combined config: tree parameters + optional seed/line overrides.
# If lv/rv seeds are not provided, base_seed is used for both.



# --- Default input paths ---
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
FILES_ORGANIZE = ROOT / "filesOrganize"
INPUTS = FILES_ORGANIZE / "cardiac_preproc" / "inputs" / "meshes"
CONFIGS = ROOT / "cardiac_preproc" / "configs"

MESH = str(INPUTS / "biv_ellipsoid.msh")
# Surface tag name mapping (used by the fractal generator).
SURFACE_TAGS = {
    "ENDO_LV": "endo_surface_lv",
    "ENDO_RV": "endo_surface_rv",
    "EPI": "epi_surface",
}

# Transmural constraint (point data field + max allowed value).
transmural_field = "uvc_transmural"
transmural_max = 0.8

# Output base filenames (written under outputs/)
lv_filename = "biv-line-lv"
rv_filename = "biv-line-rv"

# Base seed (used if lv_seed/rv_seed are not set)
base_seed = (24.4853, 86.5019, 371.475)

# Branch seed points (nearest mesh vertex is used)
lv_seed = (20.097, 106.584, 387.231)
rv_seed = (18.967, 94.4257, 393.146)

# Line end points for each branch (direction = end - seed)
lv_line_end = (37.0679,87.7706, 375.444)
rv_line_end = (48.374, 67.1098, 363.955)


# Initial direction vectors (used when no line is defined)
#lv_init_dir = (1.0, 0.0, 0.0)
#rv_init_dir = (1.0, 0.0, 0.0)

# LV parameters
lv_init_length = 40
lv_N_it = 12
lv_length = 10
lv_l_segment = 1
lv_branch_angle = 0.25
lv_repulsitivity = 0.1

# RV parameters
rv_init_length = 35
rv_N_it = 14
rv_length = 9.5
rv_l_segment = 1
rv_branch_angle = 0.2
rv_repulsitivity = 0.1

# Save intermediate frames
lv_save_frames = True
rv_save_frames = True

# Optional video rendering after fractal generation (pvpython).
VIDEO = {
    "enabled": True,
    "pvpython": "pvpython",
    "script": str(
        ROOT
        / "purkinje"
        / "fractal_3d"
        / "videoEditorPurkinje"
        / "render_purkinje_frames.py"
    ),
}




# Debug logging for branch growth
debug_grow = True
debug_dir = "debug_attempts"
debug_every = 100

# Plot debug points in 3D viewer (can be very slow for large trees, use with debug_grow)
plot_max_points = 2000
