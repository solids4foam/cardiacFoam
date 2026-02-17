"""Division schemes and color sets for LV angle-based segmentation."""

# Active scheme name
division_reference = "apical4_mid5_basal5"

# Each reference defines:
# - ring_layout: [septal_parts, freewall_parts] per ring
# - ring_id_map: explicit segment IDs per ring (len = septal_parts + freewall_parts)
# - apex_id: segment ID for apex cap (if enabled)
# - colors: palette list used by segment ID order (modulo if needed)
# - z_division_color/theta_division_color: overlay colors
division_references = {
    "aha17": {
        "ring_layout": {
            "basal": [2, 4],
            "mid": [2, 4],
            "apical": [1, 3],
        },
        "ring_id_map": {
            "basal": [1, 2, 3, 4, 5, 6],
            "mid": [7, 8, 9, 10, 11, 12],
            "apical": [13, 14, 15, 16],
        },
        "apex_id": 17,
        "use_apex_cap": True,
        "colors": [
            "purple", "cyan", "yellow", "magenta", "lime", "orange",
            "red", "green", "blue", "pink", "brown", "gray",
            "gold", "navy", "teal", "olive", "black",
        ],
        # (2) Coronary territory colors (LAD/LCx/RCA)
        "artery_colors": {
            "LAD": "#D62728",   # red
            "LCx": "#1F77B4",   # blue
            "RCA": "#2CA02C",   # green
        },

        # Segment -> artery mapping (typical right-dominant AHA17 convention)
        # NOTE: territory definitions can vary with dominance & patient anatomy.
        "segment_artery_map": {
            1: "LAD",  2: "LCx",  3: "LCx",  4: "RCA",  5: "RCA",  6: "LAD",
            7: "LAD",  8: "LCx",  9: "LCx", 10: "RCA", 11: "RCA", 12: "LAD",
            13: "LAD", 14: "LCx", 15: "RCA", 16: "LAD",
            17: "LAD",
        },

        "z_division_color": "white",
        "theta_division_color": "black",
    },
    # 16-segment style: no dedicated apex cap, full angular division in apical region.
    "aha16_no_apex": {
        "ring_layout": {
            "basal": [2, 4],
            "mid": [2, 4],
            "apical": [1, 3],
        },
        "ring_id_map": {
            "basal": [1, 2, 3, 4, 5, 6],
            "mid": [7, 8, 9, 10, 11, 12],
            "apical": [13, 14, 15, 16],
        },
        "apex_id": None,
        "use_apex_cap": False,
        "colors": [
            "purple", "cyan", "yellow", "magenta", "lime", "orange",
            "red", "green", "blue", "pink", "brown", "gray",
            "gold", "navy", "teal", "olive",
        ],
        "z_division_color": "white",
        "theta_division_color": "black",
    },
    
    # Template: apical 4 sectors, mid 5 sectors, basal 5 sectors (+ optional apex cap) V. Garcia-Bustos et al.2017
    "apical4_mid5_basal5": {
        "ring_layout": {
            "basal": [2, 3],   # total 5
            "mid": [2, 3],     # total 5
            "apical": [1, 3],  # total 4
        },
        "ring_id_map": {
            "basal": [1, 2, 3, 4, 5],
            "mid": [6, 7, 8, 9, 10],
            "apical": [11, 12, 13, 14],
        },
        "apex_id": None,
        "use_apex_cap": False,
        "colors": [
            "purple", "cyan", "yellow", "magenta", "lime",
            "orange", "red", "green", "blue", "pink",
            "brown", "gray", "gold", "navy",
        ],
        "z_division_color": "white",
        "theta_division_color": "black",
    },
}
