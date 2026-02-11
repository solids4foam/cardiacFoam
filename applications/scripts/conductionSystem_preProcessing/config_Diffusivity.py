"""
config.py

Configuration for Purkinje layer, conversion/inspection, and diffusivity.
"""

# --- Purkinje layer defaults (layer mode) ---

    #Layer thickness as fraction of transmural depth( 10% transmural depth starting in subendocardium)
PURKINJE_TRANSMURAL_MIN = 0.0
PURKINJE_TRANSMURAL_MAX = 0.1

    #Tags for LV and RV in uvc_intraventricular field (Check fields list to confirm)
LV_INTRAVENTRICULAR_TAG = -1
RV_INTRAVENTRICULAR_TAG =  1

    #Named for the cell data field tagging Purkinje elements
PURKINJE_FIELD_NAME = "purkinjeLayer"



# --- Diffusivities (ventricular tissue) in SI units: S/m ---
Ventricular_DF = 0.107  # Along fiber direction
Ventricular_DS = 0.049  # Along sheet direction
Ventricular_DN = 0.016  # Normal to sheet direction

# --- Diffusivity in purkinje Layer scale multiplier ---
PURKINJE_DIFFUSION_MULTIPLIER = 2.0

