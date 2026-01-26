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



    #Named for the cell data field tagging Purkinje elements for the slab approach
PURKINJE_FIELD_NAME = "purkinjeLayer"
