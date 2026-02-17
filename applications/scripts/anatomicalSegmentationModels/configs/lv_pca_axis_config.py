"""Plane values for lv_pca_axis plotting."""



z_apical_mid = 0.4
z_mid_basal = 0.6
z_cap_valves = 0.8

z_apical_apex = 0.2

LV_INTRAVENTRICULAR_TAG = -1
RV_INTRAVENTRICULAR_TAG = 1

# Anatomy-based circumferential segmentation inputs.
# Set groove angles directly (radians), or leave as None and provide mask fields.
ant_groove_theta = None
post_groove_theta = None

# Optional point-data field names for groove masks (nonzero == True).
ant_groove_field = True
post_groove_field = True

# Optional side-by-side reference image shown with angle-division plot.
angle_reference_image = "aha_lv17_reference.png"

# Optional fallback point on septal wall to orient grooves when RV reference is unavailable.
# Format: (x, y, z) or None.
septal_click_xyz = None
