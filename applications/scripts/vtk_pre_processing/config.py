"""
config.py

Configuration file for mesh processing pipeline.
Defines file paths, default diffusivities, and other constants.
"""

# --- File paths ---

# Folders
OUTPUT_FOLDER = "output"


INPUT_FILE = "ASCIIlegacy.vtk"
MIDDLE_FILE = f"{OUTPUT_FOLDER}/purkinjeNoDiffusion.vtk"
OUTPUT_FILE = f"{OUTPUT_FOLDER}/Banana.vtk"


# --- Default diffusivities (ventricular tissue) in SI units: S/m ---
Ventricular_DF =0.107  # Along fiber direction
Ventricular_DS = 0.049  # Along sheet direction
Ventricular_DN = 0.016  # Normal to sheet direction
